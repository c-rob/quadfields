// This file was created for V-REP release V3.4.0 

#include "libv_repExtSymDeriv.hpp"

#define CONCAT(x,y,z) x y z
#define strConCat(x,y,z)    CONCAT(x,y,z)
#define EX_TO_DOUBLE(x)	GiNaC::ex_to<numeric>(x).to_double();

using namespace GiNaC;
using namespace std;

LIBRARY vrepLib; // the V-REP library that we will dynamically load and bind


// Plugin global variables
unsigned nVars = 0;					// the number of symbolic vars (up to 4)

// dynamic properties
float mass = 0;
matrix J_inertia = {{1,0,0},{0,1,0},{0,0,1}};


// 	The flat output derivatives 
ex flatOut[4];			// numeric values in the current iteration
ex flatOut1[4];
ex flatOut2[4];
ex flatOut3[4];
ex flatOut4[4];


// State vector; saved in vrep conventions (angles and frames)
struct State {
	double x, y, z, vx, vy, vz, a, b, g, p, q, r;
};

// Symbolic equations
struct {
	ex phi;
	ex theta;
	ex psi;
	ex d_phi;
	ex d_theta;
	ex d_psi;
	matrix T;
	matrix omega;
	matrix d_omega;
	matrix u_torque;
	ex u_thrust;
} equations;

symbol Sx("x"), Sy("y"), Sz("z"), Syaw("w");	// the variables
symbol St("t"); 		// Time is implicit in all variables
matrix flatOut_D1;		// symbolic flat output derivatives vectors
matrix flatOut_D2;
matrix flatOut_D3;
matrix flatOut_D4;


// forward declaration
void updateState(Inputs &inputs, double x, double y, double z, double yaw);


/*
 Declare symbolic functions for Ginac authomatic differentiation
*/
static ex diff_symF(const ex &var, const ex &nDiff, const ex &t, unsigned diff_param);
static ex eval_symF(const ex &var, const ex &nDiff, const ex &t);
static ex evalf_symF(const ex &var, const ex &nDiff, const ex &t);

// Symbolic function for generic variables
DECLARE_FUNCTION_3P(symF)
REGISTER_FUNCTION(symF, evalf_func(evalf_symF).
		derivative_func(diff_symF).
		eval_func(eval_symF))


static ex diff_symF(const ex &var, const ex &nDiff, const ex &t, unsigned diff_param) {
	if (diff_param == 2) {
		return symF(var, nDiff+1, t);
	} else if (diff_param == 0) {
		return 1;
	} else {
		cerr << "symF Bad differentiation\n" << endl;
		return 0;
	}
}


static ex eval_symF(const ex &var, const ex &nDiff, const ex &t) {

	// Simplify symbolic equations if the vector field do not specifies z or yaw
	if (var.is_equal(Sz) && nVars < 3 && nDiff > 0) {
		return 0;
	} else if (var.is_equal(Syaw) && nVars < 4 && nDiff > 0) {
		return 0;
	} else {
		return symF(var, nDiff, t).hold();
	}
}


static ex evalf_symF(const ex &var, const ex &nDiff, const ex &t) {
	// Warning: flat outputs must be evaluated first

	unsigned index = 0;
	if (var.is_equal(Sx)) index = 0;
	else if (var.is_equal(Sy)) index = 1;
	else if (var.is_equal(Sz)) index = 2;
	else if (var.is_equal(Syaw)) index = 3;
	else {
		cerr << "Error: wrong index in symF\n" << endl;
		return 0;
	}

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
		cerr << "Error: wrong symbol as argument of symF\n" << endl;
		return 0;
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
    3,
	sim_script_arg_string,1,
	sim_script_arg_double,1,
	sim_script_arg_table | sim_script_arg_double,9
};
const int inArgs_UPDATE[]={
    4,
    sim_script_arg_double,1,
    sim_script_arg_double,1,
    sim_script_arg_double,1,
    sim_script_arg_double,1
};


void LUA_INIT_CALLBACK(SScriptCallBack* cb)
{ 
    CScriptFunctionData D;
	int ret = false;
    if (D.readDataFromStack(cb->stackID,inArgs_INIT,inArgs_INIT[0],LUA_INIT_COMMAND))
    {
		// fileName
        std::vector<CScriptFunctionDataItem>* inData=D.getInDataPtr();
		string fileName = inData->at(0).stringData[0];

		// mass
		mass = inData->at(1).doubleData[0];

		// inertia matrix: TODO: move inside initField and rotate
		for (unsigned r = 0; r < 3; ++r) {
			for (unsigned c = 0; c < 3; ++c) {
				J_inertia(r,c) = inData->at(2).doubleData[r*3+c];
			}
		}


		// call
		ret = initField(fileName);
    }
    D.pushOutData(CScriptFunctionDataItem(ret));
    D.writeDataToStack(cb->stackID);
}

 

void LUA_UPDATE_CALLBACK(SScriptCallBack* cb)
{ 
    CScriptFunctionData D;
	Inputs inputs;
    if (D.readDataFromStack(cb->stackID,inArgs_UPDATE,inArgs_UPDATE[0],LUA_UPDATE_COMMAND))
    {
        std::vector<CScriptFunctionDataItem>* inData=D.getInDataPtr();
		double x = inData->at(0).doubleData[0];
		double y = inData->at(1).doubleData[0];
		double z = inData->at(2).doubleData[0];
		double yaw = inData->at(3).doubleData[0];

		// call
		updateState(inputs, x, y, z, yaw);

    }
	// return quadrotor inputs
    D.pushOutData(CScriptFunctionDataItem(inputs.fz));
    D.pushOutData(CScriptFunctionDataItem(inputs.tx));
    D.pushOutData(CScriptFunctionDataItem(inputs.ty));
    D.pushOutData(CScriptFunctionDataItem(inputs.tz));
    D.writeDataToStack(cb->stackID);
}
// --------------------------------------------------------------------------------------


matrix vectorVrepTransform(const matrix& vec) {

	return matrix({{vec(0,0)}, {-vec(1,0)}, {-vec(2,0)}});
}


matrix rpy2matrix(const ex r, const ex p, const ex y) {
	// [r,p,y] = [phi,theta,psi]

	matrix Rz = {{cos(y), -sin(y), 0}, {sin(y), cos(y), 0}, {0, 0, 1}};
	matrix Ry = {{cos(p), 0, sin(p)}, {0, 1, 0}, {-sin(p), 0, cos(p)}};
	matrix Rx = {{1, 0, 0}, {0, cos(r), -sin(r)}, {0, sin(r), cos(r)}};

	return Rz.mul(Ry.mul(Rx));
}


matrix abg2matrix(const ex a, const ex b, const ex g) {
	// [a,b,g] = [alpha,beta,gamma]

	matrix Rz = {{cos(g), -sin(g), 0}, {sin(g), cos(g), 0}, {0, 0, 1}};
	matrix Ry = {{cos(b), 0, sin(b)}, {0, 1, 0}, {-sin(b), 0, cos(b)}};
	matrix Rx = {{1, 0, 0}, {0, cos(a), -sin(a)}, {0, sin(a), cos(a)}};

	return Rx.mul(Ry.mul(Rz));
}


matrix matrix2rpy(matrix m) {
	// [r,p,y] = [phi,theta,psi]

	ex phi = atan2(m(2,1), m(2,2));
	ex theta = atan2(-m(2,0), sqrt(m(2,2)*m(2,2) + m(2,1)*m(2,1)));
	ex psi = atan2(m(1,0), m(0,0));

	return matrix({{phi}, {theta}, {psi}});
}


matrix matrix2abg(matrix m) {
	// [a,b,g] = [alpha,beta,gamma]

	ex a = atan2(-m(1,2), m(2,2));
	ex b = atan2(m(0,2), sqrt(m(1,2)*m(1,2) + m(2,2)*m(2,2)));
	ex g = atan2(-m(0,1), m(0,0));

	return matrix({{a}, {b}, {g}});
}


// given the vector of the variables, the symbolic vector src in inputs
// computes the next derivative through dv/dt = J_v(x) * dx/dt
//  NOTE: This is different from the reference paper! They wrote dv/dt = J_v(x) * v
void genNextDerivative(const vector <symbol> vars, const matrix& src,
        const matrix& dx, matrix& dest) {

	unsigned nVars = vars.size();
	matrix jacob(nVars, nVars);
	for (unsigned r = 0; r < nVars; ++r) {
		for (unsigned c = 0; c < nVars; ++c) {
			jacob(r, c) = src[r].diff(vars.at(c));
		}
	}

	dest = jacob.mul(dx);
}


void flatOutputs2state(State &state) {

	// Endogenous transformation: state in paper, eq.8
	// state: x, y, z, vx , vy , vz , psi, theta, phi, p, q, r

	ex x = flatOut[0];
	ex y = flatOut[1];
	ex z = flatOut[2];

	ex vx = flatOut1[0];
	ex vy = flatOut1[1];
	ex vz = flatOut1[2];

	ex phi = equations.phi.evalf();
	ex theta = equations.theta.evalf();
	ex psi = flatOut[3];

	ex p = equations.omega(0,0).evalf();
	ex q = equations.omega(1,0).evalf();
	ex r = equations.omega(2,0).evalf();
	 
	
	// Transform to Vrep convention
	matrix vTemp = vectorVrepTransform(matrix({{state.x},{state.y},{state.z}}));
	state.x = EX_TO_DOUBLE(vTemp(0,0));
	state.y = EX_TO_DOUBLE(vTemp(1,0));
	state.z = EX_TO_DOUBLE(vTemp(2,0));

	vTemp = vectorVrepTransform(matrix({{state.vx},{state.vy},{state.vz}}));
	state.vx = EX_TO_DOUBLE(vTemp(0,0));
	state.vy = EX_TO_DOUBLE(vTemp(1,0));
	state.vz = EX_TO_DOUBLE(vTemp(2,0));

	matrix abg = matrix2abg(rpy2matrix(phi, theta, psi));
	state.a = EX_TO_DOUBLE(abg(0,0));
	state.b = EX_TO_DOUBLE(abg(1,0));
	state.g = EX_TO_DOUBLE(abg(2,0));

	vTemp = vectorVrepTransform(matrix({{state.p},{state.q},{state.r}}));
	state.p = EX_TO_DOUBLE(vTemp(0,0));
	state.q = EX_TO_DOUBLE(vTemp(1,0));
	state.r = EX_TO_DOUBLE(vTemp(2,0));

	// Debug
	cout << "state: " << state.x << ", " << state.y << ", " << state.z << "; "
	 	<< state.vx << ", " << state.vy << ", " << state.vz << "; "
	 	<< state.a << ", " << state.b << ", " << state.g << "; "
	 	<< state.p << ", " << state.q << ", " << state.r << endl;
}


void flatOutputs2inputs(Inputs &inputs) {

	// TODO: convert in vrep convention? inertia?

	// substitute symbolic equations
	ex u_torque_f = equations.u_torque.evalf();
	
	inputs.tx = EX_TO_DOUBLE(ex_to<matrix>(u_torque_f)(0,0));
	inputs.ty = EX_TO_DOUBLE(ex_to<matrix>(u_torque_f)(1,0));
	inputs.tz = EX_TO_DOUBLE(ex_to<matrix>(u_torque_f)(2,0));
	inputs.fz = EX_TO_DOUBLE(equations.u_thrust.evalf());

}


void genSymbolicEquations(void) {

	// Euler angles and first derivative
	ex ba = -cos(symF(Syaw,0,St)) * symF(Sx,2,St) - sin(symF(Syaw,0,St)) * symF(Sy,2,St);
	ex bb = -symF(Sz,2,St) + 9.8;
	ex bc = -sin(symF(Syaw,0,St)) * symF(Sx,2,St) + cos(symF(Syaw,0,St)) * symF(Sy,2,St);

	equations.phi = atan2(bc, sqrt(ba*ba + bb*bb));
	equations.theta = atan2(ba, bb);
	equations.psi = symF(Syaw,0,St);

	equations.d_phi = equations.phi.diff(St);
	equations.d_theta = equations.theta.diff(St);
	equations.d_psi = equations.psi.diff(St);

	// Euler rates rpy to angular velocity
	equations.T = {
		{ cos(equations.phi) * cos(equations.theta), -sin(equations.phi), 0},
		{ cos(equations.theta) * sin(equations.phi), cos(equations.phi), 0},
		{ -sin(equations.theta), 0, 1}
	};
	ex temp_omega = equations.T * matrix({{equations.d_phi},
			{equations.d_theta}, {equations.d_psi}});
	ex temp_d_omega = temp_omega.diff(St);

	equations.omega = ex_to<matrix>(temp_omega.evalm());
	equations.d_omega = ex_to<matrix>(temp_d_omega.evalm());


	// omega = [0, −r, q; r, 0, −p; −q, p, 0]
	matrix skewOmega = {{0, -equations.omega(2,0), equations.omega(1,0)},
						{equations.omega(2,0), 0, -equations.omega(0,0)},
						{-equations.omega(1,0), equations.omega(0,0), 0}};

	//equations.u_torque = J_inertia * equations.d_omega + skew(equations.omega)
	//	* J_inertia * equations.omega;
	ex temp_u_torque = J_inertia * equations.d_omega + skewOmega * J_inertia * equations.omega;
	equations.u_torque = ex_to<matrix>(temp_u_torque.evalm());
	
	// equations.u_thrust = m_mass * norm(flatOut_D2[0:2] - 9.8 * [0;0;1])
	matrix xyz_D2 = {{symF(Sx,2,St)},{symF(Sy,2,St)},{symF(Sz,2,St)}};
	//ex flatOut_D2_sube = (xyz_D2 - matrix({{0},{0},{9.8}}));
	ex flatOut_D2_sube = (xyz_D2);
	matrix flatOut_D2_subm = ex_to<matrix>(flatOut_D2_sube.evalm());
	ex innerProd = flatOut_D2_subm.transpose() * flatOut_D2_subm;
	matrix innerProdM = ex_to<matrix>(innerProd.evalm());
	equations.u_thrust = mass * sqrt(innerProdM(0,0));
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
	
	// Save the first flat output derivative d(sigma)/dt=V(x)
	flatOut_D1 = vectFieldSym;

	// Compute next derivatives
	genNextDerivative(vars, flatOut_D1, flatOut_D1, flatOut_D2);
	genNextDerivative(vars, flatOut_D2, flatOut_D1, flatOut_D3);
	genNextDerivative(vars, flatOut_D3, flatOut_D1, flatOut_D4);

	// Save equations to globals
	genSymbolicEquations();

	return true;
}


// The registered vrep function for evaluating the inputs
void updateState(Inputs &inputs, double x, double y, double z, double yaw) {
	
	// Pass from the v-rep axis convention to reference paper conv. (z downwards)
	matrix vTemp = vectorVrepTransform(matrix({{x}, {y}, {z}}));
	x = EX_TO_DOUBLE(vTemp(0,0));
	y = EX_TO_DOUBLE(vTemp(1,0));
	z = EX_TO_DOUBLE(vTemp(2,0));
	// NOTE: yaw is ignored! (assumed 0)
	// TODO: check and fix axis conversions

	// Evaluate the D4 vectors numerically
	exmap symMap;
	symMap[Sx] = x;
	symMap[Sy] = y;
	symMap[Sz] = z;
	symMap[Syaw] = yaw;

	// fill the globals flatOutputs derivatives
	for (unsigned i = 0; i < nVars; ++i) {
		flatOut1[i] = flatOut_D1(i,0).subs(symMap).evalf();
		flatOut2[i] = flatOut_D2(i,0).subs(symMap).evalf();
		flatOut3[i] = flatOut_D3(i,0).subs(symMap).evalf();
		flatOut4[i] = flatOut_D4(i,0).subs(symMap).evalf();
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
	flatOutputs2state(state);
	flatOutputs2inputs(inputs);

	// DEBUG
	//cout << "vecField: " << flatOut1[0] << ", " << flatOut1[1] << ", " <<
	//	flatOut1[2] << ", " << flatOut1[3] << endl;
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
			strConCat("number ok = ",LUA_INIT_COMMAND,"(string filePath, number mass, table9 inertia_matrix)"),LUA_INIT_CALLBACK);

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

	// CLN setting of 0 exception
	cln::cl_inhibit_floating_point_underflow = true;

	initField("/home/roberto/Desktop/Erob/V-REP/symsplugin/vector-field.txt");
	/*
	mass = 1;

	Inputs inp;
	updateState(inp, 2, 3, 4, 0);
	*/


}

