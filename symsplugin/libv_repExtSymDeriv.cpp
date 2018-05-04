// This file was created for V-REP release V3.4.0 

#include "libv_repExtSymDeriv.hpp"


// #define DEBUG
// #define DEBUG_PRINT
// #define DEBUG_PRINT_FLAT_OUTPUTS


#define CONCAT(x,y,z) x y z
#define strConCat(x,y,z)    CONCAT(x,y,z)
#define EX_TO_DOUBLE(x)	GiNaC::ex_to<numeric>(x).to_double()

using namespace GiNaC;
using namespace std;


LIBRARY vrepLib; // the V-REP library that we will dynamically load and bind

// GinaC settings
bool cln::cl_inhibit_floating_point_underflow = true; // no underflow exception

/***
 * Globals
 ***/
int quadcopterH = -1;
unsigned nVars = 0;					// The number of lines in the vector field file

// dynamic properties
float mass = 0;
matrix J_inertia = {{1,0,0},{0,1,0},{0,0,1}};


/***
 * Numeric quantities in the current iteration
 ***/
ex flatOut[4];			// measured flat outputs
ex flatOut1[4];			// computed derivatives:
ex flatOut2[4];
ex flatOut3[4];
ex flatOut4[4];

// State vector; saved in vrep conventions (angles and frames)
struct State {
	double x, y, z, vx, vy, vz, a, b, g, p, q, r;
};


/***
 * Symbolic equation, computed at initialization
 ***/
symbol Sx("x"), Sy("y"), Sz("z"), Syaw("w");	// the variables (flat outputs)
symbol St("t"); 		// Time is implicit in all variables

struct {
	ex phi;				// Rotation about x
	ex theta;			// Rotation about y
	ex psi;				// Rotation about z
	ex d_phi;			// ^ deriv
	ex d_theta; 		// ^ deriv
	ex d_psi;   		// ^ deriv
	matrix omega;		// Angular vel in body frame
	matrix d_omega;		// ^ deriv
	matrix u_torque;	// Control input: torque in x,y,z
	ex u_thrust;		// Control input: thrust

	matrix R;
	matrix d_R;
} equations;

matrix flatOut_D1;		// vectors of symbolic flat output derivatives
matrix flatOut_D2;
matrix flatOut_D3;
matrix flatOut_D4;



// Debug symbols to delete
ex intVel = matrix(3, 1);
ex intAbg = matrix(3, 1);
ex intRpy = matrix(3, 1);
ex abgLast = matrix(3, 1);


/*
 * Code starts
 */

// forward declaration
void updateState(Inputs &inputs, State &state, double x, double y, double z,
        double a, double b, double g);

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

	// TODO: is this an error? reference at 0 doesn't mean flatOut to 0
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
	// NOTE: flat outputs must be evaluated before any .evalf()!

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
	sim_script_arg_table | sim_script_arg_double,9,
};
const int inArgs_UPDATE[]={
    2,
	sim_script_arg_table | sim_script_arg_double,3,
	sim_script_arg_table | sim_script_arg_double,3,
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

		// inertia matrix
		for (unsigned r = 0; r < 3; ++r) {
			for (unsigned c = 0; c < 3; ++c) {
				J_inertia(r,c) = inData->at(2).doubleData[r*3+c];
			}
		}

		// call
		ret = initField(fileName, true);
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
		double y = inData->at(0).doubleData[1];
		double z = inData->at(0).doubleData[2];
		double a = inData->at(1).doubleData[0];
		double b = inData->at(1).doubleData[1];
		double g = inData->at(1).doubleData[2];

		// call
	    State state;
		updateState(inputs, state, x, y, z, a, b, g);

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


matrix rpy2matrix(matrix rpy) {
	// [r,p,y] = [phi,theta,psi]
	
	ex r = rpy(0,0);
	ex p = rpy(1,0);
	ex y = rpy(2,0);

	matrix Rz = {{cos(y), -sin(y), 0}, {sin(y), cos(y), 0}, {0, 0, 1}};
	matrix Ry = {{cos(p), 0, sin(p)}, {0, 1, 0}, {-sin(p), 0, cos(p)}};
	matrix Rx = {{1, 0, 0}, {0, cos(r), -sin(r)}, {0, sin(r), cos(r)}};

	return Rz.mul(Ry.mul(Rx));
}


matrix abg2matrix(matrix abg) {
	// [a,b,g] = [alpha,beta,gamma]

	ex a = abg(0,0);
	ex b = abg(1,0);
	ex g = abg(2,0);

	matrix Rz = {{cos(g), -sin(g), 0}, {sin(g), cos(g), 0}, {0, 0, 1}};
	matrix Ry = {{cos(b), 0, sin(b)}, {0, 1, 0}, {-sin(b), 0, cos(b)}};
	matrix Rx = {{1, 0, 0}, {0, cos(a), -sin(a)}, {0, sin(a), cos(a)}};

	// NOTE: Should I take the inverted axis into account?
	matrix Ryi = Ry.transpose();
	matrix Rzi = Rz.transpose();

	return Rx.mul(Ryi.mul(Rzi));
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

	// NOTE: Opposite b,g: should I take the inverted axis into account?
	ex a = atan2(-m(1,2), m(2,2));
	ex b = atan2(-m(0,2), sqrt(m(1,2)*m(1,2) + m(2,2)*m(2,2)));
	ex g = atan2(m(0,1), m(0,0));

	return matrix({{a}, {b}, {g}});
}


matrix omega2rpyRate(matrix angVel, matrix rpy) {

	matrix TInv = {
		{ 1, (sin(rpy(0,0))*sin(rpy(1,0)))/cos(rpy(1,0)), (cos(rpy(0,0))*sin(rpy(1,0)))/cos(rpy(1,0))},
		{ 0, cos(rpy(0,0)), -sin(rpy(0,0))},
		{ 0, sin(rpy(0,0))/cos(rpy(1,0)), cos(rpy(0,0))/cos(rpy(1,0))}};

	return TInv.mul(angVel);
}


matrix rpyRate2omega(matrix rpyRate, matrix rpy) {

	matrix T = {
		{ 1, 0, -sin(rpy(1,0))},
		{ 0, cos(rpy(0,0)), cos(rpy(1,0))*sin(rpy(0,0))},
		{ 0, -sin(rpy(0,0)), cos(rpy(0,0))*cos(rpy(1,0))}};

	return T.mul(rpyRate);
}


matrix omega2abgRate(matrix angVelVrep, matrix abg) {

	// NOTE: changing reference is enough?
	matrix angVel = vectorVrepTransform(angVelVrep);

	matrix TInv = {
		{ cos(abg(2,0))/cos(abg(1,0)), -sin(abg(2,0))/cos(abg(1,0)), 0},
		{ sin(abg(2,0)), cos(abg(2,0)), 0},
		{ -(cos(abg(2,0))*sin(abg(1,0)))/cos(abg(1,0)), (sin(abg(1,0))*sin(abg(2,0)))/cos(abg(1,0)), 1}};

	return TInv.mul(angVel);
}


matrix abgRate2omega(matrix abgRate, matrix abg) {

	matrix T = {
		{ cos(abg(1,0))*cos(abg(2,0)), sin(abg(2,0)), 0},
		{ -cos(abg(1,0))*sin(abg(2,0)), cos(abg(2,0)), 0},
		{ sin(abg(1,0)), 0, 1}};

	matrix omegaVrep = T.mul(abgRate);

	// NOTE: changing reference is enough?
	return vectorVrepTransform(omegaVrep);
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
	matrix vTemp = vectorVrepTransform(matrix({{x},{y},{z}}));
	state.x = EX_TO_DOUBLE(vTemp(0,0));
	state.y = EX_TO_DOUBLE(vTemp(1,0));
	state.z = EX_TO_DOUBLE(vTemp(2,0));

	vTemp = vectorVrepTransform(matrix({{vx},{vy},{vz}}));
	state.vx = EX_TO_DOUBLE(vTemp(0,0));
	state.vy = EX_TO_DOUBLE(vTemp(1,0));
	state.vz = EX_TO_DOUBLE(vTemp(2,0));

        // This must be compared to the dummy object in vrep scene (paper axis convention)
	matrix abg = matrix2abg(rpy2matrix(matrix({{phi}, {theta}, {psi}})));
	state.a = EX_TO_DOUBLE(ex_to<matrix>(abg)(0,0));
	state.b = EX_TO_DOUBLE(ex_to<matrix>(abg)(1,0));
	state.g = EX_TO_DOUBLE(ex_to<matrix>(abg)(2,0));

	vTemp = vectorVrepTransform(matrix({{p},{q},{r}}));
	state.p = EX_TO_DOUBLE(ex_to<matrix>(vTemp)(0,0));
	state.q = EX_TO_DOUBLE(ex_to<matrix>(vTemp)(1,0));
	state.r = EX_TO_DOUBLE(ex_to<matrix>(vTemp)(2,0));

}


void flatOutputs2inputs(Inputs &inputs) {

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

	// Euler rates rpy to angular velocity (remeber: the result is omega in local frame)
	equations.omega = rpyRate2omega(matrix({{equations.d_phi},{equations.d_theta},{equations.d_psi}}),
			matrix({{equations.phi},{equations.theta},{equations.psi}}));

	ex temp_d_omega = equations.omega.diff(St);
	equations.d_omega = ex_to<matrix>(temp_d_omega.evalm());


	// omega = [0, −r, q; r, 0, −p; −q, p, 0]  (in local frame too)
	matrix skewOmega = {{0, -equations.omega(2,0), equations.omega(1,0)},
						{equations.omega(2,0), 0, -equations.omega(0,0)},
						{-equations.omega(1,0), equations.omega(0,0), 0}};

	//equations.u_torque = J_inertia * equations.d_omega + skew(equations.omega)
	//	* J_inertia * equations.omega;
	ex temp_u_torque = J_inertia * equations.d_omega + skewOmega * J_inertia * equations.omega;
	equations.u_torque = ex_to<matrix>(temp_u_torque.evalm());

	// equations.u_thrust = m_mass * norm(flatOut_D2[0:2] - 9.8 * [0;0;1])
	matrix xyz_D2 = {{symF(Sx,2,St)},{symF(Sy,2,St)},{symF(Sz,2,St)}};
	ex flatOut_D2_sube = (xyz_D2 - matrix({{0},{0},{9.8}}));
	//ex flatOut_D2_sube = (xyz_D2); // NOTE: with or without gravity compensation?
	matrix flatOut_D2_subm = ex_to<matrix>(flatOut_D2_sube.evalm());
	ex innerProd = flatOut_D2_subm.transpose() * flatOut_D2_subm;
	matrix innerProdM = ex_to<matrix>(innerProd.evalm());
	equations.u_thrust = mass * sqrt(innerProdM(0,0));

	
	// Additional equations useful for debugging
	equations.R = rpy2matrix(matrix({{equations.phi},{equations.theta},{equations.psi}}));
	equations.d_R = ex_to<matrix>(equations.R.diff(St));
}


void setVrepInitialState(void) {


    // get initial position of the quadcopter shape
    quadcopterH = simGetObjectHandle("Quadricopter");
	float initPos[3];
	simGetObjectPosition(quadcopterH, -1, initPos);
	
	// Compute a configuration in the field
    Inputs inputs;
    State state;
    updateState(inputs, state, initPos[0], initPos[1], initPos[2], 0, 0, 0);
        // NOTE: arg 8 is the 4-th flat output, set to 0 here. 6-7 args are not needed
	
	// Set other properties
	const float abg[] = {(float)state.a, (float)state.b, (float)state.g};
	simSetObjectOrientation(quadcopterH, -1, abg);

	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_x, (float)state.vx);
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_y, (float)state.vy);
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_z, (float)state.vz);

	matrix abgD1 = omega2abgRate(matrix({{state.p}, {state.q}, {state.r}}), matrix({{state.a}, {state.b}, {state.g}}));
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_a, EX_TO_DOUBLE(abgD1(0,0)));
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_b, EX_TO_DOUBLE(abgD1(1,0)));
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_g, EX_TO_DOUBLE(abgD1(2,0)));

#ifdef DEBUG_PRINT
	cout << "initPos " << initPos[0] << ", " << initPos[1] << ", " << initPos[2] << endl;
	cout << "initVel " << state.vx << ", " << state.vy << ", " << state.vz << endl;
	cout << "orient " << abg[0] << ", " << abg[1] << ", " << abg[2] << endl;
	cout << "Dorient " << abgD1 << endl;
#endif

#ifdef DEBUG
	intVel = vectorVrepTransform(matrix({{state.vx}, {state.vy}, {state.vz}}));
	intAbg = matrix({{abg[0]},{abg[1]},{abg[2]}});
	intRpy = matrix(3,1); // assuming start from horizontal
#endif
}


int initField(string fieldFilePath, bool vrepCaller) {

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
	//vars.erase(vars.begin()+nVars, vars.end());
		// NOTE: if this is commented, differentiation is on all 4 vars (usually nothing changes)

	// Fill a symbolic matrix
	matrix vectFieldSym(4, 1);
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

    // Assigns initial config in vrep scene to match the vector field
	if (vrepCaller) {
		setVrepInitialState();
	}

	return true;
}


void debugging(Inputs &inputs, State &state, double x, double y, double z,
        double a, double b, double g) {

	// state
	matrix pos = {{x}, {y}, {z}};
	matrix velV = {{state.vx},{state.vy},{state.vz}};
	matrix vel = vectorVrepTransform(velV);
	matrix acc = {{flatOut2[0]},{flatOut2[1]},{flatOut2[2]}};
	matrix jerk = {{flatOut3[0]},{flatOut3[1]},{flatOut3[2]}};
	matrix snap = {{flatOut4[0]},{flatOut4[1]},{flatOut4[2]}};
	matrix abg = {{state.a}, {state.b}, {state.g}};
	matrix pqr = {{state.p}, {state.q}, {state.r}};
	matrix vrepAbg = {{a}, {b}, {g}};

	// integrate velocities
	ex newPos = pos + 0.01 * vel;
	matrix newPosM = ex_to<matrix>(newPos.evalm());
	matrix newPosVM = vectorVrepTransform(newPosM);
	float newPosVF[3];
	newPosVF[0] = EX_TO_DOUBLE(newPosVM(0,0));
	newPosVF[1] = EX_TO_DOUBLE(newPosVM(1,0));
	newPosVF[2] = EX_TO_DOUBLE(newPosVM(2,0));

	// integrate deriv rpy
	matrix dRpy = {{(equations.d_phi.evalf())},
		{(equations.d_theta.evalf())},
		{(equations.d_psi.evalf())}};
	cout << "dRpy " << dRpy << endl;
	intRpy = (intRpy + 0.01 * dRpy).evalm();
	cout << "intRpy " << intRpy << endl;
	matrix intAbg = matrix2abg(rpy2matrix(ex_to<matrix>(intRpy)));
	

	float abgF[3];
	abgF[0] = EX_TO_DOUBLE(abg(0,0));
	abgF[1] = EX_TO_DOUBLE(abg(1,0));
	abgF[2] = EX_TO_DOUBLE(abg(2,0));
	simSetObjectOrientation(quadcopterH, -1, abgF);
	simSetObjectPosition(quadcopterH, -1, newPosVF);

	cout << "int.abg:   " << intAbg << endl;;
	cout << "state.abg: " << abg << endl << endl;

	// debug: off motors
	inputs.fz = 0;
	inputs.tx = 0;
	inputs.ty = 0;
	inputs.tz = 0;
}


// The registered vrep function for evaluating the inputs
void updateState(Inputs &inputs, State &state, double x, double y, double z,
        double a, double b, double g) {

	// Pass from the v-rep axis convention to reference paper conv. (z downwards)
	matrix vTemp = vectorVrepTransform(matrix({{x}, {y}, {z}}));
	x = EX_TO_DOUBLE(vTemp(0,0));
	y = EX_TO_DOUBLE(vTemp(1,0));
	z = EX_TO_DOUBLE(vTemp(2,0));

    matrix rpy = matrix2rpy(abg2matrix(matrix({{a}, {b}, {g}})));
	double yaw = EX_TO_DOUBLE(rpy(2,0));

	// Evaluate the D4 vectors numerically
	exmap symMap;
	symMap[Sx] = x;
	symMap[Sy] = y;
	symMap[Sz] = z;
	symMap[Syaw] = yaw;

	// fill the globals flatOutputs derivatives
	for (unsigned i = 0; i < 4; ++i) {
		flatOut1[i] = flatOut_D1(i,0).subs(symMap).evalf();
		flatOut2[i] = flatOut_D2(i,0).subs(symMap).evalf();
		flatOut3[i] = flatOut_D3(i,0).subs(symMap).evalf();
		flatOut4[i] = flatOut_D4(i,0).subs(symMap).evalf();
	}

	// save to global
	flatOut[0] = x;
	flatOut[1] = y;
	flatOut[2] = z;
	flatOut[3] = yaw;

#ifdef DEBUG_PRINT_FLAT_OUTPUTS
	// Deb_print
	cout << "flatOut: " << flatOut[0] << ", " << flatOut[1] << ", " <<
		flatOut[2] << ", " << flatOut[3] << endl;
	cout << "flatOut1: " << flatOut1[0] << ", " << flatOut1[1] << ", " <<
		flatOut1[2] << ", " << flatOut1[3] << endl;
	cout << "flatOut2: " << flatOut2[0] << ", " << flatOut2[1] << ", " <<
		flatOut2[2] << ", " << flatOut2[3] << endl;
	cout << "flatOut3: " << flatOut3[0] << ", " << flatOut3[1] << ", " <<
		flatOut3[2] << ", " << flatOut3[3] << endl;
	cout << "flatOut4: " << flatOut4[0] << ", " << flatOut4[1] << ", " <<
		flatOut4[2] << ", " << flatOut4[3] << endl;
#endif


	// Get the state of the quadrotor
	flatOutputs2state(state);
	flatOutputs2inputs(inputs);

#ifdef DEBUG
	debugging(inputs, state, x, y, z, a, b, g);
#endif
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
	
	// init properties
	mass = 0.87;
	J_inertia.add(ex_to<matrix>(diag_matrix({0.006,0.006,0.012})));

	initField("./vector-field.txt", false);

	// set a fictitious pose
	Inputs inputs;
	State state;
	updateState(inputs, state, 1,0,1, 0,0,1);

	// Print initial state
	matrix abg = {{state.a}, {state.b}, {state.g}};
	matrix pqr = {{state.p}, {state.q}, {state.r}};
	cout << "Init state:\n";
	cout << "pos:  " << state.x << ", " << state.y << ", " << state.z << endl;
	cout << "vel:  " << state.vx << ", " << state.vy << ", " << state.vz << endl;
	cout << "abg:  " << state.a << ", " << state.b << ", " << state.g << endl;
	cout << "pqr:  " << state.p << ", " << state.q << ", " << state.r << endl << endl;

	// phi,th,psi,R,d_phi,d_th,d_psi equations are ok

	// Evaluate all equations
	matrix rpy = {{equations.phi.evalf()},{equations.theta.evalf()},{equations.psi.evalf()}};
	matrix d_rpy = {{equations.d_phi.evalf()},{equations.d_theta.evalf()},{equations.d_psi.evalf()}};
	ex omega = equations.omega.evalf();
	ex d_omega = equations.d_omega.evalf();
	ex R = equations.R.evalf();
	ex d_R = equations.d_R.evalf();
	ex Omega = (equations.R.transpose().evalf() * equations.d_R.evalf()).evalm();

	cout << "\n\nNumeric\n";
	cout << "rpy:   " << rpy << endl;
	cout << "d_rpy: " << d_rpy << endl;
	cout << "omega: " << omega << endl;
	cout << "d_omega: " << d_omega << endl;
	cout << "R:     " << R << endl;
	cout << "d_R:   " << d_R << endl;
	cout << "Omega:     " << Omega << endl;


	// DEBUG !!
	// 		+ Check equations.omega comparing with d_R * R.transpose()
	// 		+ check we can differentiate angular velocities directly (theory)
	// 		+ compare d_omega with its expression with diff rotation matrices
	// 		+ if necessary compute omega and its derivatives with diff rot
	//
}

