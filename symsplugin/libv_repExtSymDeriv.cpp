// This file was created for V-REP release V3.4.0 

#include "libv_repExtSymDeriv.hpp"

#define CONCAT(x,y,z) x y z
#define strConCat(x,y,z)    CONCAT(x,y,z)

using namespace GiNaC;
using namespace std;

LIBRARY vrepLib; // the V-REP library that we will dynamically load and bind


// Plugin global variables
unsigned nVars = 0;					// the number of symbolic vars (up to 4)
symbol Sx("x"), Sy("y"), Sz("z"), Syaw("w");	// the variables
symbol St("t"); 		// Time is implicit in all variables
matrix flatOut_D1;		// symbolic flat output derivatives vectors
matrix flatOut_D2;
matrix flatOut_D3;
matrix flatOut_D4;

// 	The flat output derivatives 
ex flatOut[4];			// numeric values in the current iteration
ex flatOut1[4];
ex flatOut2[4];
ex flatOut3[4];
ex flatOut4[4];


// state vector
struct State {
	ex x, y, z, vx, vy, vz, psi, theta, phi, p, q, r;
};
// input vector: thrust + torques
struct Inputs {
	ex fz, tx, ty, tz;
};

// Symbolic equations
struct {
	ex phi;
	ex theta;
	ex psi;
	ex d_phi;
	ex d_theta;
	ex d_psi;
} equations;



/*
 Declare symbolic functions for Ginac authomatic differentiation
*/
static ex diff_symF(const ex &var, const ex &nDiff, const ex &t, unsigned diff_param);
static ex evalf_symF(const ex &var, const ex &nDiff, const ex &t);

// Symbolic function for generic variables
DECLARE_FUNCTION_3P(symF)
REGISTER_FUNCTION(symF, evalf_func(evalf_symF). derivative_func(diff_symF))


static ex diff_symF(const ex &var, const ex &nDiff, const ex &t, unsigned diff_param) {
	if (diff_param == 2) {
		return symF(var, nDiff+1, t);
	} else if (diff_param == 0) {
		return 1;
	} else {
		throw "symF Bad differentiation\n";
	}
}


static ex evalf_symF(const ex &var, const ex &nDiff, const ex &t) {
	// Warning: flat outputs must be evaluated first

	unsigned index = 0;
	if (var.is_equal(Sx)) index = 0;
	else if (var.is_equal(Sy)) index = 1;
	else if (var.is_equal(Sz)) index = 2;
	else if (var.is_equal(Syaw)) index = 3;
	else throw "Error: wrong index in symF\n";

	ex ret;
	numeric nDiffN = ex_to<numeric>(nDiff);
	if (nDiffN == 0) {
		ret = flatOut[index];
	} else if (nDiffN == 1) {
		ret = flatOut1[index];
	} else if (nDiffN == 2) {
		ret = flatOut2[index];
	} else if (nDiffN == 3) {
		ret = flatOut3[index];
	} else if (nDiffN == 4) {
		ret = flatOut4[index];
	} else {
		throw "Error: wrong symbol as argument of symF\n";
	}

	return ret;
}



// --------------------------------------------------------------------------------------
// simExtSymDeriv_init
// --------------------------------------------------------------------------------------
#define LUA_INIT_COMMAND "simExtSymDeriv_init"
 
// --------------------------------------------------------------------------------------
// simExtSymDeriv_update
// --------------------------------------------------------------------------------------
#define LUA_UPDATE_COMMAND "simExtSymDeriv_update"

const int inArgs_INIT[]={
    1,
	sim_script_arg_string,1
};
const int inArgs_UPDATE[]={
    4,
    sim_script_arg_float,1,
    sim_script_arg_float,1,
    sim_script_arg_float,1,
    sim_script_arg_float,1
};


void LUA_INIT_CALLBACK(SScriptCallBack* cb)
{ 
    CScriptFunctionData D;
	int ret = false;
    if (D.readDataFromStack(cb->stackID,inArgs_INIT,inArgs_INIT[0],LUA_INIT_COMMAND))
    {
        std::vector<CScriptFunctionDataItem>* inData=D.getInDataPtr();
		string fileName = inData->at(0).stringData[0];

		// call
		ret = initField(fileName);
    }
    D.pushOutData(CScriptFunctionDataItem(ret));
    D.writeDataToStack(cb->stackID);
}

 

void LUA_UPDATE_CALLBACK(SScriptCallBack* cb)
{ 
    CScriptFunctionData D;
    if (D.readDataFromStack(cb->stackID,inArgs_UPDATE,inArgs_UPDATE[0],LUA_UPDATE_COMMAND))
    {
        std::vector<CScriptFunctionDataItem>* inData=D.getInDataPtr();
		float x = inData->at(0).floatData[0];
		float y = inData->at(1).floatData[0];
		float z = inData->at(2).floatData[0];
		float yaw = inData->at(3).floatData[0];

		// call
		updateState(x, y, z, yaw);
    }
	// return quadrotor inputs // TODO: correct values
    D.pushOutData(CScriptFunctionDataItem(3.3));
    D.pushOutData(CScriptFunctionDataItem(3.4));
    D.pushOutData(CScriptFunctionDataItem(3.5));
    D.pushOutData(CScriptFunctionDataItem(3.6));
    D.writeDataToStack(cb->stackID);
}
// --------------------------------------------------------------------------------------


// given the vector of the variables, the symbolic vector src in inputs
// computes the next derivative through dv = J_v * v
void genNextDerivative(const vector <symbol> vars, const matrix& src, matrix& dest) {

	unsigned nVars = vars.size();
	matrix jacob(nVars, nVars);
	for (unsigned r = 0; r < nVars; ++r) {
		for (unsigned c = 0; c < nVars; ++c) {
			jacob(r, c) = src[r].diff(vars.at(c));
		}
	}

	dest = jacob.mul(src);
}


void flatOutputs2state(State &state, const ex flatOut[], const ex flatOut1[],
		const ex flatOut2[], const ex flatOut3[]) {

	// Endogenous transformation: state in paper, eq.8
	// state: x, y, z, vx , vy , vz , psi, theta, phi, p, q, r

	state.x = flatOut[0];
	state.y = flatOut[1];
	state.z = flatOut[2];
	state.vx = flatOut1[0];
	state.vy = flatOut1[1];
	state.vz = flatOut1[2];

	ex ba = -cos(flatOut[3]) * flatOut2[0] - sin(flatOut[3]) * flatOut2[1];
	ex bb = -flatOut2[2] + 9.8;
	ex bc = -sin(flatOut[3]) * flatOut2[0] + cos(flatOut[3]) * flatOut2[1];

	state.phi = atan2(bc, sqrt(ba*ba + bb*bb));
	state.theta = atan2(ba, bb);
	state.psi = flatOut[3];

	// Euler rates rpy to angular velocity
	matrix T = {
		{ cos(state.phi) * cos(state.theta), -sin(state.phi), 0},
		{ cos(state.theta) * sin(state.phi), cos(state.phi), 0},
		{ -sin(state.theta), 0, 1}
	};

	// computing derivatives for [p q r] computations:
	
	// d_ba = sin(s4(t))*diff(s4(t), t)*diff(s1(t), t, t) - sin(s4(t))*diff(s2(t), t, t, t) + 
	// - cos(s4(t))*diff(s4(t), t)*diff(s2(t), t, t) - cos(s4(t))*diff(s1(t), t, t, t)
	ex d_ba = sin(flatOut[3])*flatOut1[3]*flatOut2[0] - sin(flatOut[3])*flatOut3[1]
	- cos(flatOut[3])*flatOut1[3]*flatOut2[1] - cos(flatOut[3])*flatOut3[0];

	ex d_bb = -flatOut3[2];

	// d_bc = cos(s4(t))*diff(s2(t), t, t, t) - sin(s4(t))*diff(s1(t), t, t, t) +
	// - cos(s4(t))*diff(s4(t), t)*diff(s1(t), t, t) - sin(s4(t))*diff(s4(t), t)*diff(s2(t), t, t)
	ex d_bc = cos(flatOut[3])*flatOut3[1] - sin(flatOut[3])*flatOut3[0]
	- cos(flatOut[3])*flatOut1[3]*flatOut2[0] - sin(flatOut[3])*flatOut1[3]*flatOut2[1];

	// d_theta = d_Atan(ba(t)/bb(t))
	ex d_theta = 1 / (1 + (ba / bb)*(ba / bb)) * (d_ba * bb - ba * d_bb) / (bb * bb);

	// d_phi = (diff(bc_(t), t)/(ba_(t)^2 + bb_(t)^2)^(1/2) - (bc_(t)*(ba_(t)*diff(ba_(t), t)
	// + bb_(t)*diff(bb_(t), t)))/(ba_(t)^2 + bb_(t)^2)^(3/2))/(bc_(t)^2/(ba_(t)^2 + bb_(t)^2) + 1)
	ex d_phi = (d_bc/pow(ba*ba + bb*bb, 1/2) - (bc*(ba*d_ba
		+ bb*d_bb))/pow(ba*ba + bb*bb, 3/2))/(bc*bc/(ba*ba + bb*bb) + 1);
 
	ex d_psi = flatOut1[3];

	matrix angVel = T.mul({{d_phi}, {d_theta}, {d_psi}});
	state.p = angVel(0,0);
	state.q = angVel(1,0);
	state.r = angVel(2,0);

	// debug
	//cout << "flatOut, flatOut1:" << endl;
	//cout << flatOut[0] << ", "<< flatOut[1] << ", "<< flatOut[2] << ", "<< flatOut[3]<< endl;
	//cout << flatOut1[0] << ", "<< flatOut1[1] << ", "<< flatOut1[2] << ", "<< flatOut1[3]<< endl;

	//cout << "state:" << endl;
	//cout << state.x << ", "<< state.y << ", "<< state.z << endl;
	//cout << state.vx << ", "<< state.vy << ", "<< state.vz << endl;
	//cout << state.phi << ", "<< state.theta << ", "<< state.psi << endl;
	//cout << state.p << ", "<< state.q << ", "<< state.r << endl;

}


void flatOutputs2inputs(Inputs &inputs, const ex flatOut[], const ex flatOut1[],
		const ex flatOut2[]) {

	// computing algular acceleration via:
	// 	d(omega)/dt = d(T(phi,theta) * [d(phi); d(theta); d(psi)])/dt



}


void genSymbolicEquations(void) {

	ex ba = -cos(symF(Syaw,0,St)) * symF(Sx,2,St) - sin(symF(Syaw,0,St)) * symF(Sy,2,St);
	ex bb = -symF(Sz,2,St) + 9.8;
	ex bc = -sin(symF(Syaw,0,St)) * symF(Sx,2,St) + cos(symF(Syaw,0,St)) * symF(Sy,2,St);

	equations.phi = atan2(bc, sqrt(ba*ba + bb*bb));
	equations.theta = atan2(ba, bb);
	equations.psi = symF(Syaw,0,St);

	equations.d_phi = equations.phi.diff(St);
	equations.d_theta = equations.theta.diff(St);
	equations.d_psi = equations.psi.diff(St);

}


int initField(string fieldFilePath) {

	string line;
	ifstream vectFile;
	vector<string> vectFieldStr;

	// File open
	try {
		vectFile.open(fieldFilePath, ifstream::in);
	} catch (ifstream::failure e) {
		cout << e.what();
		return false;
	}

	// Prepare the GiNaC parser
	symtab table;
	vector <symbol> vars = {Sx, Sy, Sz, Syaw};
	table["x"] = Sx;
	table["y"] = Sy;
	table["z"] = Sz;
	table["w"] = Syaw;
	parser reader(table);

	// Get the first lines in the file as a vector
	nVars = 0;
	while (getline(vectFile, line)) {
		vectFieldStr.push_back(line);
		++nVars;
	}
	vectFile.close();
	vars.erase(vars.begin()+nVars, vars.end());

	// Fill a symbolic matrix
	matrix vectFieldSym(nVars, 1);
	for (unsigned i = 0; i < nVars; ++i) {
		ex e = reader(vectFieldStr[i]);
		vectFieldSym.set(i, 0, e);
	}

	// Pass from the v-rep axis convention to reference paper conv. (z downwards)
	if (nVars >= 2) {
		vectFieldSym.set(1, 0, -vectFieldSym(1, 0));
	}
	if (nVars >= 3) {
		vectFieldSym.set(2, 0, -vectFieldSym(2, 0));
	}

	// Save the first flat output derivative d(sigma)/dt=V(x)
	flatOut_D1 = vectFieldSym;

	// Compute next derivatives
	genNextDerivative(vars, flatOut_D1, flatOut_D2);
	genNextDerivative(vars, flatOut_D2, flatOut_D3);
	genNextDerivative(vars, flatOut_D3, flatOut_D4);

	// Save equations to globals
	genSymbolicEquations();

	return true;
}


// The registered vrep function for evaluating the inputs
void updateState(float x, float y, float z, float yaw) {
	
	// Pass from the v-rep axis convention to reference paper conv. (z downwards)
	y = -y;
	z = -z;

	// Evaluate the D4 vectors numerically
	exmap symMap;
	symMap[Sx] = x;
	symMap[Sy] = y;
	symMap[Sz] = z;
	symMap[Syaw] = yaw;

	// fill the globals flatOutputs derivatives
	for (unsigned i = 0; i < nVars; ++i) {
		flatOut1[i] = flatOut_D1(i,0).subs(symMap);
		flatOut2[i] = flatOut_D2(i,0).subs(symMap);
		flatOut3[i] = flatOut_D3(i,0).subs(symMap);
		flatOut4[i] = flatOut_D4(i,0).subs(symMap);
		//flatOut1[i] = ex_to<numeric>(flatOut_D1(i,0).subs(symMap)).to_double();
		//flatOut2[i] = ex_to<numeric>(flatOut_D2(i,0).subs(symMap)).to_double();
		//flatOut3[i] = ex_to<numeric>(flatOut_D3(i,0).subs(symMap)).to_double();
		//flatOut4[i] = ex_to<numeric>(flatOut_D4(i,0).subs(symMap)).to_double();
	}
	for (unsigned i = nVars; i < 4; ++i) {
		flatOut1[i] = 0;
		flatOut2[i] = 0;
		flatOut3[i] = 0;
		flatOut4[i] = 0;
	}
	// save to global
	flatOut[0] = x;
	flatOut[1] = y;
	flatOut[2] = z;
	flatOut[3] = yaw;

	// Get the state of the quadrotor
	State state;
	flatOutputs2state(state, flatOut, flatOut1, flatOut2, flatOut3);

}



// This is the plugin start routine (called just once, just after the plugin was loaded):
VREP_DLLEXPORT unsigned char v_repStart(void* reservedPointer,int reservedInt)
{
    // Dynamically load and bind V-REP functions:
    // 1. Figure out this plugin's directory:
    char curDirAndFile[1024];
#ifdef _WIN32
    #ifdef QT_COMPIL
        _getcwd(curDirAndFile, sizeof(curDirAndFile));
    #else
        GetModuleFileName(NULL,curDirAndFile,1023);
        PathRemoveFileSpec(curDirAndFile);
    #endif
#else
    getcwd(curDirAndFile, sizeof(curDirAndFile));
#endif

    std::string currentDirAndPath(curDirAndFile);
    // 2. Append the V-REP library's name:
    std::string temp(currentDirAndPath);
#ifdef _WIN32
    temp+="\\v_rep.dll";
#elif defined (__linux)
    temp+="/libv_rep.so";
#elif defined (__APPLE__)
    temp+="/libv_rep.dylib";
#endif /* __linux || __APPLE__ */
    // 3. Load the V-REP library:
    vrepLib=loadVrepLibrary(temp.c_str());
    if (vrepLib==NULL)
    {
        std::cout << "Error, could not find or correctly load the V-REP library. Cannot start 'PluginSkeleton' plugin.\n";
        return(0); // Means error, V-REP will unload this plugin
    }
    if (getVrepProcAddresses(vrepLib)==0)
    {
        std::cout << "Error, could not find all required functions in the V-REP library. Cannot start 'PluginSkeleton' plugin.\n";
        unloadVrepLibrary(vrepLib);
        return(0); // Means error, V-REP will unload this plugin
    }

    // Check the version of V-REP:
    int vrepVer;
    simGetIntegerParameter(sim_intparam_program_version,&vrepVer);
    if (vrepVer<30200) // if V-REP version is smaller than 3.02.00
    {
        std::cout << "Sorry, your V-REP copy is somewhat old. Cannot start 'PluginSkeleton' plugin.\n";
        unloadVrepLibrary(vrepLib);
        return(0); // Means error, V-REP will unload this plugin
    }

    // Register the lua commands
    simRegisterScriptCallbackFunction(strConCat(LUA_INIT_COMMAND,"@","SymDeriv"),
			strConCat("number ok = ",LUA_INIT_COMMAND,"(string filePath)"),LUA_INIT_CALLBACK);

    simRegisterScriptCallbackFunction(strConCat(LUA_UPDATE_COMMAND,"@","SymDeriv"),
			strConCat("",LUA_UPDATE_COMMAND,"(number x, number y, number z, number yaw)"),LUA_UPDATE_CALLBACK);


    return(PLUGIN_VERSION); // initialization went fine, we return the version number of this plugin (can be queried with simGetModuleName)
}

// This is the plugin end routine (called just once, when V-REP is ending, i.e. releasing this plugin):
VREP_DLLEXPORT void v_repEnd()
{
    // Here you could handle various clean-up tasks

    unloadVrepLibrary(vrepLib); // release the library
}

// This is the plugin messaging routine (i.e. V-REP calls this function very often, with various messages):
VREP_DLLEXPORT void* v_repMessage(int message,int* auxiliaryData,void* customData,int* replyData)
{ // This is called quite often. Just watch out for messages/events you want to handle
    // Keep following 5 lines at the beginning and unchanged:
    static bool refreshDlgFlag=true;
    int errorModeSaved;
    simGetIntegerParameter(sim_intparam_error_report_mode,&errorModeSaved);
    simSetIntegerParameter(sim_intparam_error_report_mode,sim_api_errormessage_ignore);
    void* retVal=NULL;

    // Here we can intercept many messages from V-REP (actually callbacks). Only the most important messages are listed here.
    // For a complete list of messages that you can intercept/react with, search for "sim_message_eventcallback"-type constants
    // in the V-REP user manual.

    if (message==sim_message_eventcallback_refreshdialogs)
        refreshDlgFlag=true; // V-REP dialogs were refreshed. Maybe a good idea to refresh this plugin's dialog too

    if (message==sim_message_eventcallback_menuitemselected)
    { // A custom menu bar entry was selected..
        // here you could make a plugin's main dialog visible/invisible
    }

    if (message==sim_message_eventcallback_instancepass)
    {   // This message is sent each time the scene was rendered (well, shortly after) (very often)
        // It is important to always correctly react to events in V-REP. This message is the most convenient way to do so:

        int flags=auxiliaryData[0];
        bool sceneContentChanged=((flags&(1+2+4+8+16+32+64+256))!=0); // object erased, created, model or scene loaded, und/redo called, instance switched, or object scaled since last sim_message_eventcallback_instancepass message 
        bool instanceSwitched=((flags&64)!=0);

        if (instanceSwitched)
        {
            // React to an instance switch here!!
        }

        if (sceneContentChanged)
        { // we actualize plugin objects for changes in the scene

            //...

            refreshDlgFlag=true; // always a good idea to trigger a refresh of this plugin's dialog here
        }
    }

    if (message==sim_message_eventcallback_mainscriptabouttobecalled)
    { // The main script is about to be run (only called while a simulation is running (and not paused!))
        
    }

    if (message==sim_message_eventcallback_simulationabouttostart)
    { // Simulation is about to start

    }

    if (message==sim_message_eventcallback_simulationended)
    { // Simulation just ended

    }

    if (message==sim_message_eventcallback_moduleopen)
    { // A script called simOpenModule (by default the main script). Is only called during simulation.
        //if ( (customData==NULL)||(_stricmp("PluginSkeleton",(char*)customData)==0) ) // is the command also meant for this plugin?
        //{
        //    // we arrive here only at the beginning of a simulation
        //}
    }

    if (message==sim_message_eventcallback_modulehandle)
    { // A script called simHandleModule (by default the main script). Is only called during simulation.
        //if ( (customData==NULL)||(_stricmp("PluginSkeleton",(char*)customData)==0) ) // is the command also meant for this plugin?
        //{
        //    // we arrive here only while a simulation is running
        //}
    }

    if (message==sim_message_eventcallback_moduleclose)
    { // A script called simCloseModule (by default the main script). Is only called during simulation.
        //if ( (customData==NULL)||(_stricmp("PluginSkeleton",(char*)customData)==0) ) // is the command also meant for this plugin?
        //{
        //    // we arrive here only at the end of a simulation
        //}
    }

    if (message==sim_message_eventcallback_instanceswitch)
    { // We switched to a different scene. Such a switch can only happen while simulation is not running

    }

    if (message==sim_message_eventcallback_broadcast)
    { // Here we have a plugin that is broadcasting data (the broadcaster will also receive this data!)

    }

    if (message==sim_message_eventcallback_scenesave)
    { // The scene is about to be saved. If required do some processing here (e.g. add custom scene data to be serialized with the scene)

    }

    // You can add many more messages to handle here

    if ((message==sim_message_eventcallback_guipass)&&refreshDlgFlag)
    { // handle refresh of the plugin's dialogs
        // ...
        refreshDlgFlag=false;
    }

    // Keep following unchanged:
    simSetIntegerParameter(sim_intparam_error_report_mode,errorModeSaved); // restore previous settings
    return(retVal);
}



int main() {
	initField("/home/roberto/Desktop/Erob/V-REP/symsplugin/vector-field.txt");

	updateState(1,-2,3,0);
	cout << equations.d_phi.evalf() << endl;

}

