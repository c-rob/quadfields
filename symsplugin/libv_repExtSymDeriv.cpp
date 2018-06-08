// This file was created for V-REP release V3.4.0 

#include "libv_repExtSymDeriv.hpp"


#define DEBUG
#define DEBUG_PRINT_INIT
// #define DEBUG_PRINT_FLAT_OUTPUTS
#define DEBUG_PRINT_INPUTS


#define CONCAT(x,y,z) x y z
#define strConCat(x,y,z)	CONCAT(x,y,z)
#define EX_TO_DOUBLE(x)	GiNaC::ex_to<numeric>(x).to_double()
#define GINAC_3VEC(x) matrix({{x[0]}, {x[1]}, {x[2]}})

using namespace GiNaC;
using namespace std;


LIBRARY vrepLib; // the V-REP library that we will dynamically load and bind

// GinaC settings
bool cln::cl_inhibit_floating_point_underflow = true; // no underflow exception

/***
 * Globals
 ***/
int quadcopterH = -1;				// vrep handle
unsigned nVars = 0;					// The number of lines in the vector field file

// dynamic properties
float mass = 0;
matrix J_inertia = {{1,0,0},{0,1,0},{0,0,1}};
const float GRAVITY_G = 9.80655;


/***
 * Numeric quantities in the current iteration
 ***/
ex flatOut[4];			// measured flat outputs
ex flatOut1[4];			// computed derivatives:
ex flatOut2[4];
ex flatOut3[4];
ex flatOut4[4];

// State vector; saved in vrep conventions (angles and frames)
//		p,q,r is the angular velocity in global vrep frame
struct State {
	double x, y, z, vx, vy, vz, a, b, g, p, q, r;
};

// Input vecotor struct defined in header


/***
 * Symbolic equation, computed at initialization
 ***/
symbol Sx("x"), Sy("y"), Sz("z"), Syaw("w");	// the variables (flat outputs)
symbol St("t");			// Time is implicit in all variables

// All these quantities are in rpy paper convention
struct {
	ex phi;				// Rotation about x
	ex theta;			// Rotation about y
	ex psi;				// Rotation about z
	ex d_phi;			// ^ deriv
	ex d_theta;			// ^ deriv
	ex d_psi;			// ^ deriv
	matrix omega;		// Angular vel in body frame
	matrix d_omega;		// ^ deriv
	matrix u_torque;	// Control input: torque in x,y,z
	ex u_thrust;		// Control input: thrust

	matrix R;
	matrix d_R;
	matrix dd_R;
} equations;

matrix flatOut_D1;		// vectors of symbolic flat output derivatives
matrix flatOut_D2;
matrix flatOut_D3;
matrix flatOut_D4;


// Debug variables for integration
const float dt = 0.005;
TinyIntegrator linVelInt("(paper) Linear veocity", dt),	// 3x1 vector
			   linPosInt("(paper) Position", dt),		// 3x1 vector
			   d_RInt("(paper) d_R", dt),				// 3x3 matrix
			   RInt("(paper) R", dt);					// 3x3 matrix
simFloat oldVrepMatrix[12];


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

	// NOTE: the first two ifs should never happen; not well tested.
	//	How to derive the last flat outputs if the field is unspecified for them?
	
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

	// NOTE: this must be coherent with vectorVrepR()

	return matrix({{vec(0,0)}, {-vec(1,0)}, {-vec(2,0)}});
}


matrix vectorVrepR() {

	// NOTE: this must be coherent with vectorVrepR()
	
	return matrix({{1,0,0},{0,-1,0},{0,0,-1}});
}


matrix matrixVrepTransform(const matrix& mat) {

	// This is a transformation that returns the matrix corresponding to
	// the same rotation matrix, in the other axis convention
	//	(equal to its inverse)
	
	matrix Rpv = vectorVrepR();		// NOTE: assuming this is equal to its inverse
	ex R = Rpv * mat * Rpv;
	return ex_to<matrix>(R.evalm());
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
	//
	// NOTE: all from, to vrep axis: classic rotation

	ex a = abg(0,0);
	ex b = abg(1,0);
	ex g = abg(2,0);

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

	// NOTE: all from, to vrep axis: classic rotation

	ex a = atan2(-m(1,2), m(2,2));
	ex b = atan2(m(0,2), sqrt(m(0,1)*m(0,1) + m(0,0)*m(0,0)));
	ex g = atan2(-m(0,1), m(0,0));

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

// This is the standart T to get omega in fixed frame, using the upper version
//	to get omega in body frame
//		{ cos(rpy(2,0))*cos(rpy(1,0)), -sin(rpy(2,0)), 0 },
//		{ sin(rpy(2,0))*cos(rpy(1,0)), cos(rpy(2,0)), 0 },
//		{ -sin(rpy(1,0)), 0, 1 }};

	return T.mul(rpyRate);
}


matrix omegaVrep2abgRate(matrix angVelVrep, matrix abg) {

	// NOTE: omega in vrep global frame

	ex a = abg(0,0);
	ex b = abg(1,0);

	matrix TInv = {
		{ 1, (sin(a)*sin(b))/cos(b), -(cos(a)*sin(b))/cos(b)},
		{ 0,				 cos(a),				  sin(a)},
		{ 0,		 -sin(a)/cos(b),		   cos(a)/cos(b)}
	};
/*
 *
 *	  // This is the matrix for omega in vrep body frame
 *
 *	  ex b = abg(1,0);
 *	  ex g = abg(2,0);
 *
 *	  matrix TInv = {
 *		  {			  cos(g)/cos(b),		 -sin(g)/cos(b), 0},
 *		  {					 sin(g),				 cos(g), 0},
 *		  { -(cos(g)*sin(b))/cos(b), (sin(b)*sin(g))/cos(b), 1}
 *	  };
 *
 */

	return TInv.mul(angVelVrep);

}


matrix abgRate2omegaVrep(matrix abgRate, matrix abg) {
	
	// NOTE: omega in vrep global frame

	ex a = abg(0,0);
	ex b = abg(1,0);

	matrix T = {
		{ 1,	  0,		 sin(b) },
		{ 0, cos(a), -cos(b)*sin(a) },
		{ 0, sin(a),  cos(a)*cos(b) }
	};

/*
 *	  // This is the matrix for omega in vrep body frame
 *
 *	  ex b = abg(1,0);
 *	  ex g = abg(2,0);
 *
 *	  matrix T = {
 *		  { cos(b)*cos(g), sin(g), 0},
 *		  {-cos(b)*sin(g), cos(g), 0},
 *		  { sin(b)		 ,		0, 1}
 *	  };
 *
 */
	return T.mul(abgRate);

}


// given the vector of the variables, the symbolic vector src in inputs
// computes the next derivative through dv/dt = J_v(x) * dx/dt
//	NOTE: This is different from the reference paper! They wrote dv/dt = J_v(x) * v
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

	// Compute numeric values for all equations regarding state quantities

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

	ex omegaGlobal = equations.R * equations.omega;
	matrix omegaGlobalM = ex_to<matrix>(omegaGlobal.evalm().evalf());

	// Transform to Vrep convention
	matrix vTemp = vectorVrepTransform(matrix({{x},{y},{z}}));
	state.x = EX_TO_DOUBLE(vTemp(0,0));
	state.y = EX_TO_DOUBLE(vTemp(1,0));
	state.z = EX_TO_DOUBLE(vTemp(2,0));

	vTemp = vectorVrepTransform(matrix({{vx},{vy},{vz}}));
	state.vx = EX_TO_DOUBLE(vTemp(0,0));
	state.vy = EX_TO_DOUBLE(vTemp(1,0));
	state.vz = EX_TO_DOUBLE(vTemp(2,0));

	matrix Rpaper = rpy2matrix(matrix({{phi}, {theta}, {psi}}));
	matrix Rvrep = matrixVrepTransform(Rpaper);
	matrix abg = matrix2abg(Rvrep);
	state.a = EX_TO_DOUBLE(abg(0,0));
	state.b = EX_TO_DOUBLE(abg(1,0));
	state.g = EX_TO_DOUBLE(abg(2,0));

	matrix vrepOmegaGlob = vectorVrepTransform(omegaGlobalM);
	state.p = EX_TO_DOUBLE(ex_to<matrix>(vrepOmegaGlob)(0,0));
	state.q = EX_TO_DOUBLE(ex_to<matrix>(vrepOmegaGlob)(1,0));
	state.r = EX_TO_DOUBLE(ex_to<matrix>(vrepOmegaGlob)(2,0));

}


void flatOutputs2inputs(Inputs &inputs) {

	// Compute numeric values for all input equations
	
	matrix u_torqueF = ex_to<matrix>(equations.u_torque.evalf());

	// body -> fixed -> fixedVrep -> bodyVrep
	matrix R = ex_to<matrix>(equations.R.evalf());
	ex u_torqueGlobal = R * u_torqueF;
	ex u_torqueGlobalVrepE = vectorVrepTransform(ex_to<matrix>(u_torqueGlobal.evalm()));
	matrix u_torqueVrepF = ex_to<matrix>(u_torqueGlobalVrepE.evalm());


	inputs.tx = EX_TO_DOUBLE(u_torqueVrepF(0,0));
	inputs.ty = EX_TO_DOUBLE(u_torqueVrepF(1,0));
	inputs.tz = EX_TO_DOUBLE(u_torqueVrepF(2,0));


	ex u_thrustF = equations.u_thrust.evalf();
	u_thrustF = (u_thrustF > 0) ? u_thrustF : 0;		// can't provide negative thrust
	inputs.fz = EX_TO_DOUBLE(u_thrustF);
}


void genSymbolicEquations(void) {

	// State equations first:
	ex ba = -cos(symF(Syaw,0,St)) * symF(Sx,2,St) - sin(symF(Syaw,0,St)) * symF(Sy,2,St);
	ex bb = -symF(Sz,2,St) + GRAVITY_G;
	ex bc = -sin(symF(Syaw,0,St)) * symF(Sx,2,St) + cos(symF(Syaw,0,St)) * symF(Sy,2,St);

	equations.phi = atan2(bc, sqrt(ba*ba + bb*bb));
	equations.theta = atan2(ba, bb);
	equations.psi = symF(Syaw,0,St);

	equations.d_phi = equations.phi.diff(St);
	equations.d_theta = equations.theta.diff(St);
	equations.d_psi = equations.psi.diff(St);
			// NOTE: ^ we shouldn't diff this but get from omega (ok for now:
			// small angles assumption)

	// Euler rates rpy to angular velocity (remeber: the result is omega in local frame)
	equations.omega = rpyRate2omega(matrix({{equations.d_phi},{equations.d_theta},{equations.d_psi}}),
			matrix({{equations.phi},{equations.theta},{equations.psi}}));

	ex temp_d_omega = equations.omega.diff(St);
	equations.d_omega = ex_to<matrix>(temp_d_omega.evalm());


	// omega = [0, −r, q; r, 0, −p; −q, p, 0]  (in local frame too)
	matrix skewOmega = {{0, -equations.omega(2,0), equations.omega(1,0)},
						{equations.omega(2,0), 0, -equations.omega(0,0)},
						{-equations.omega(1,0), equations.omega(0,0), 0}};

	// Inputs: torque
	ex temp_u_torque = J_inertia * equations.d_omega + skewOmega * J_inertia * equations.omega;
	equations.u_torque = ex_to<matrix>(temp_u_torque.evalm());

	// Inputs: thrust
		// equations.u_thrust = m_mass * norm(flatOut_D	2[0:2] - GRAVITY_G * [0;0;1])
	matrix e3 = {{0},{0},{1}};
	matrix xyzD2 = {{symF(Sx,2,St)},{symF(Sy,2,St)},{symF(Sz,2,St)}};
	ex thrustAccVec = (xyzD2 - GRAVITY_G * e3);			// NOTE: with gravity compensation?
	matrix thrustAccVecM = ex_to<matrix>(thrustAccVec.evalm());
	ex thrustAccNorm = thrustAccVecM.transpose() * thrustAccVecM;		// norm of the acceleration vector
	matrix tempMat = ex_to<matrix>(thrustAccNorm.evalm());
	equations.u_thrust = mass * sqrt(tempMat(0,0));		// thrust absolute value


	// Other useful, but unnecessary, equations
	equations.R = rpy2matrix(matrix({{equations.phi},{equations.theta},{equations.psi}}));
	equations.d_R = ex_to<matrix>(equations.R.diff(St));
	equations.dd_R = ex_to<matrix>(equations.d_R.diff(St));
	

	// These are equivalent expressions for quantities already computed
	/*
	matrix Omega = R.transpose().mul(d_R);

	ex d_OmegaEx = d_R.transpose() * d_R + R.transpose() * dd_R;
	matrix d_Omega = ex_to<matrix>(d_OmegaEx.evalm());

	matrix xyzD2 = {{symF(Sx,2,St)},{symF(Sy,2,St)},{symF(Sz,2,St)}};
	ex thrustVec = equations.R.transpose() * mass * (GRAVITY_G * e3 - xyz_D2);
	matrix fxVecM = ex_to<matrix>(thrustVec.evalm());
	equations.u_thrust = fxVecM(2,0);				// NOTE: negative values are set to 0 in flatOutputs2inputs()
	*/
}


void setVrepInitialState(void) {

	// TODO: make the name of the shape parametric

	// Get the initial pose of the quadcopter shape in the vrep scene
	quadcopterH = simGetObjectHandle("Quadricopter");
	float initVrepPos[3];
	float initVrepAbg[3];
	simGetObjectOrientation(quadcopterH, -1, initVrepAbg);
	simGetObjectPosition(quadcopterH, -1, initVrepPos);

	// Compute a configuration in the field
	Inputs inputs;
	State state;
	updateState(inputs, state, initVrepPos[0], initVrepPos[1], initVrepPos[2], initVrepAbg[0], initVrepAbg[1], initVrepAbg[2]);
		// NOTE: arg 8 is the 4-th flat output. 6-7 args can be different from state.a,state.b
	

	// This section computes R_t1 from the continuous derivative of R_t0
	// Then, it extracts abg from the delta transformation as needed by the 
	// sim_shapefloatparam_init_velocity_{a,b,g} parameters
	matrix R_t0 = ex_to<matrix>(equations.R.evalf());
	matrix d_R_t0 = ex_to<matrix>(equations.d_R.evalf());
	ex R_t1 = R_t0 + d_R_t0 * dt;
	ex R_t0t1 = R_t0.transpose() * R_t1;		// delta R
	matrix mR_t0t1 = ex_to<matrix>(R_t0t1.evalm());
	matrix vrep_R_t0t1 = matrixVrepTransform(mR_t0t1);	// applied here or to R_t{1,0} it's the same
	matrix abg_t0t1 = matrix2abg(vrep_R_t0t1);
	matrix abg_rate = abg_t0t1.mul_scalar(1.0/dt);

	// Set the vrep state
	const float abg[] = {(float)state.a, (float)state.b, (float)state.g};
	simSetObjectOrientation(quadcopterH, -1, abg);

	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_x, (float)state.vx);
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_y, (float)state.vy);
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_z, (float)state.vz);

	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_a, EX_TO_DOUBLE(abg_rate(0,0)));
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_b, EX_TO_DOUBLE(abg_rate(1,0)));
	simSetObjectFloatParameter(quadcopterH, sim_shapefloatparam_init_velocity_g, EX_TO_DOUBLE(abg_rate(2,0)));

	simResetDynamicObject(quadcopterH);

#ifdef DEBUG_PRINT_INIT
	cout << "Initital pose get:\n";
	cout << "InitVrepPos: " << initVrepPos[0] << ", " << initVrepPos[1] << ", " << initVrepPos[2] << endl;
	cout << "InitVrepAbg: " << initVrepAbg[0] << ", " << initVrepAbg[1] << ", " << initVrepAbg[2] << endl;
	cout << "Inital state set:" << endl;
	cout << "initPos " << initVrepPos[0] << ", " << initVrepPos[1] << ", " << initVrepPos[2] << endl;
	cout << "initVel " << state.vx << ", " << state.vy << ", " << state.vz << endl;
	cout << "abg " << abg[0] << ", " << abg[1] << ", " << abg[2] << endl;
	cout << "Dabg set " << abg_rate << endl;
	cout << "Other\n";
	cout << "abgD1 " << omegaVrep2abgRate(matrix({{state.p}, {state.q}, {state.r}}), matrix({{state.a}, {state.b}, {state.g}})) << endl;
	cout << "pqr " << matrix({{state.p}, {state.q}, {state.r}}) << endl << endl;
#endif

#ifdef DEBUG
	// Integrating position and rpy
	// Set integrators' initial states here
	linVelInt.setInitialState(matrix{{flatOut1[0]},{flatOut1[1]},{flatOut1[2]}});
	linPosInt.setInitialState(matrix{{flatOut[0]},{flatOut[1]},{flatOut[2]}});
	RInt.setInitialState(equations.R.evalf());		// R
	matrix omegaF = ex_to<matrix>(equations.omega.evalf());
	matrix Omega = {{0, -omegaF(2,0), omegaF(1,0)},
					{omegaF(2,0), 0, -omegaF(0,0)},
					{-omegaF(1,0), omegaF(0,0), 0}};
	d_RInt.setInitialState((equations.R.evalf() * Omega).evalm());	// d_R
	// DEBUG: checking against integration if starting still
		cout << "initInt: \n" << linVelInt << endl << d_RInt << endl << linPosInt << endl << RInt << endl << endl;
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


void debugging(Inputs& inputs, State& state) {

	// Run if initialized
	if (! RInt.isInitialized()) {
		return;
	}

	// Inputs
	ex u_thrust = (equations.u_thrust.evalf());
	matrix u_torque = ex_to<matrix>(equations.u_torque.evalf());


	// >> Begin integration: testing flat outputs to inputs to field integration
	//
	// Prepare whats needed
	matrix e3 = {{0},{0},{1}};
	matrix R = RInt.getMat();
	matrix Omega = R.transpose().mul(d_RInt.getMat());
	matrix omega = {{Omega(2,1)},{-Omega(2,0)},{Omega(1,0)}};

	// Back to accelerations
	ex accel = (GRAVITY_G * e3 - R * (u_thrust * e3 / mass)).evalm();
	ex angAcc = (J_inertia.inverse() * (u_torque - Omega * J_inertia * omega)).evalm();
	
	matrix angAccM = ex_to<matrix>(angAcc);
	matrix AngAcc = {{0, -angAccM(2,0), angAccM(1,0)},
					{angAccM(2,0), 0, -angAccM(0,0)},
					{-angAccM(1,0), angAccM(0,0), 0}};

	// Update state
	linPosInt.update(linVelInt.get());
	linVelInt.update(accel);

	d_RInt.update(R * AngAcc - R * d_RInt.getMat().transpose() * d_RInt.getMat());
	RInt.update(d_RInt.get());

	// >> End integration
	

	// >> Compare estimated and real quantities
	
	// Estimated linear velocity
	matrix vrepLinVel = matrix({{state.vx},{state.vy},{state.vz}});

	// Measured
	float vrepLinVelF[3];
	float vrepAbgVelF[3];
	simGetVelocity(quadcopterH, vrepLinVelF, vrepAbgVelF);
		// ATTENTION: at the first iteration, simGetVelocity returns just what
		// was previously set with sim_shapefloatparam_init_velocity_{a,b,g}
	
	// Computing my own angular measure
	simFloat vrepMatrix[12];
	simFloat vrepAngVelAxis[3];
	simFloat vrepAngle;
	simGetObjectMatrix(quadcopterH, -1, vrepMatrix);
	simGetRotationAxis(oldVrepMatrix, vrepMatrix, vrepAngVelAxis, &vrepAngle);
	float vrepAngVelScalar = vrepAngle/dt;
	matrix vrepAngVel = GINAC_3VEC(vrepAngVelAxis).mul_scalar(vrepAngVelScalar);
	matrix angVel = vectorVrepTransform(vrepAngVel);
	for (int i = 0; i < 12; ++i) { oldVrepMatrix[i] = vrepMatrix[i]; }

	matrix vrepOmega = vectorVrepTransform(ex_to<matrix>(equations.omega.evalf()));
	matrix vrepOmega2 = ex_to<matrix>((vrepOmega.transpose() * vrepOmega).evalm());
	ex vrepOmegaNorm = sqrt(vrepOmega2(0,0));

	// cout << "vrep est  linvel:	" << vrepLinVel << endl;
	// cout << "vrep real linvel:	" << GINAC_3VEC(vrepLinVelF) << endl;

	cout << "vrep est  axis:   " << (vrepOmega/vrepOmegaNorm).evalm() << endl;
	cout << "vrep real axis:   " << GINAC_3VEC(vrepAngVelAxis) << endl;
	cout << "vrep est angvel:  " << vrepOmegaNorm << endl;
	cout << "vrep real angvel: " << vrepAngVelScalar << endl;
	cout << endl;


	// do not set, just printing now
/*
 *    // Vrep convention
 *    matrix vrepLinPos = vectorVrepTransform(linPosInt.getMat());
 *    float vrepLinPosF[3];
 *    vrepLinPosF[0] = EX_TO_DOUBLE(vrepLinPos(0,0));
 *    vrepLinPosF[1] = EX_TO_DOUBLE(vrepLinPos(1,0));
 *    vrepLinPosF[2] = EX_TO_DOUBLE(vrepLinPos(2,0));
 *
 *    matrix vrepAbgPos = matrix2abg(matrixVrepTransform(RInt.getMat()));
 *    float vrepAbgPosF[3];
 *    vrepAbgPosF[0] = EX_TO_DOUBLE(vrepAbgPos(0,0));
 *    vrepAbgPosF[1] = EX_TO_DOUBLE(vrepAbgPos(1,0));
 *    vrepAbgPosF[2] = EX_TO_DOUBLE(vrepAbgPos(2,0));
 *
 *	  simSetObjectPosition(quadcopterH, -1, vrepLinPosF);
 *	  simSetObjectOrientation(quadcopterH, -1, vrepAbgPosF);
 *
 *	  // debug: off motors
 *	  inputs.fz = 0;
 *	  inputs.tx = 0;
 *	  inputs.ty = 0;
 *	  inputs.tz = 0;
 */
}


// The registered vrep function for evaluating the inputs
void updateState(Inputs &inputs, State &state, double x, double y, double z,
		double a, double b, double g) {

	// Pass from the v-rep axis convention to reference paper conv. (z downwards)
	matrix vTemp = vectorVrepTransform(matrix({{x}, {y}, {z}}));
	x = EX_TO_DOUBLE(vTemp(0,0));
	y = EX_TO_DOUBLE(vTemp(1,0));
	z = EX_TO_DOUBLE(vTemp(2,0));

	matrix Rvrep = abg2matrix(matrix({{a}, {b}, {g}}));
	matrix Rpaper = matrixVrepTransform(Rvrep);
	matrix rpy = matrix2rpy(Rpaper);
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

	// Get the state of the quadrotor
	flatOutputs2state(state);
	flatOutputs2inputs(inputs);

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
		flatOut4[2] << ", " << flatOut4[3] << endl << endl;
#endif

#ifdef DEBUG
	debugging(inputs, state);
#endif

#ifdef DEBUG_PRINT_INPUTS
	cout << "Inputs [fz, tx, ty, tz]: [" << inputs.fz << ", " << inputs.tx <<
		", " << inputs.ty << ", " << inputs.tz << "]\n\n";
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
	{	// This message is sent each time the scene was rendered (well, shortly after) (very often)
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
		//	  // we arrive here only at the beginning of a simulation
		//}
	}

	if (message==sim_message_eventcallback_modulehandle)
	{ // A script called simHandleModule (by default the main script). Is only called during simulation.
		//if ( (customData==NULL)||(_stricmp("PluginSkeleton",(char*)customData)==0) ) // is the command also meant for this plugin?
		//{
		//	  // we arrive here only while a simulation is running
		//}
	}

	if (message==sim_message_eventcallback_moduleclose)
	{ // A script called simCloseModule (by default the main script). Is only called during simulation.
		//if ( (customData==NULL)||(_stricmp("PluginSkeleton",(char*)customData)==0) ) // is the command also meant for this plugin?
		//{
		//	  // we arrive here only at the end of a simulation
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
	
	// Set the same Vrep dynamic properties
	mass = 0.87;
	J_inertia.set(0,0, 0.006);
	J_inertia.set(1,1, 0.006);
	J_inertia.set(2,2, 0.011);
	cout << "Mass: " << mass << endl;
	cout << "Inertia: " << J_inertia << endl;

	initField("./vector-field.txt", false);

	// set a fictitious pose
	float x = 1;
	float y = 0;
	float z = 1;
	float a = 0;
	float b = 0;
	// float g = 3.14159265359/2;
	float g = 0;
	Inputs inputs;
	State state;
	updateState(inputs, state, x, y, z, a, b, g);

	// Print flat outputs
	cout << "Flat outputs\n";
	cout << "flatOut: " << flatOut[0] << ", " << flatOut[1] << ", " <<
		flatOut[2] << ", " << flatOut[3] << endl;
	cout << "flatOut1: " << flatOut1[0] << ", " << flatOut1[1] << ", " <<
		flatOut1[2] << ", " << flatOut1[3] << endl;
	cout << "flatOut2: " << flatOut2[0] << ", " << flatOut2[1] << ", " <<
		flatOut2[2] << ", " << flatOut2[3] << endl;
	cout << "flatOut3: " << flatOut3[0] << ", " << flatOut3[1] << ", " <<
		flatOut3[2] << ", " << flatOut3[3] << endl;
	cout << "flatOut4: " << flatOut4[0] << ", " << flatOut4[1] << ", " <<
		flatOut4[2] << ", " << flatOut4[3] << "\n\n";

	// Evaluate all equations: these are well tested
	matrix rpy = {{equations.phi.evalf()},{equations.theta.evalf()},{equations.psi.evalf()}};
	matrix d_rpy = {{equations.d_phi.evalf()},{equations.d_theta.evalf()},{equations.d_psi.evalf()}};
	matrix omega = ex_to<matrix>(equations.omega.evalf());
	matrix d_omega = ex_to<matrix>(equations.d_omega.evalf());
	matrix R = ex_to<matrix>(equations.R.evalf());
	matrix d_R = ex_to<matrix>(equations.d_R.evalf());
	matrix dd_R = ex_to<matrix>(equations.dd_R.evalf());
	matrix Omega = R.transpose().mul(d_R);
	ex d_OmegaEx = d_R.transpose() * d_R + R.transpose() * dd_R;
	matrix d_Omega = ex_to<matrix>(d_OmegaEx.evalm());
	ex u_thrust = (equations.u_thrust.evalf());
	matrix u_torque = ex_to<matrix>(equations.u_torque.evalf());
	
	cout << "Paper equations\n";
	cout << "rpy:	" << rpy << endl;
	cout << "d_rpy: " << d_rpy << endl;
	cout << "omega: " << omega << endl;
	cout << "omegaGlob " << (R.mul(omega)) << endl;
	cout << "Omega: " << Omega.evalf() << endl;
	cout << "d_omega: " << d_omega << endl;
	cout << "d_omegaGlob " << (R.mul(d_omega)) << endl;
	cout << "R:		" << R << endl;
	cout << "d_R:	" << d_R << endl;
	cout << "dd_R:	 " << dd_R << endl;
	cout << "Omega:		" << Omega << endl;
	cout << "d_Omega:	" << d_Omega << endl;
	cout << "u_thrust:	 " << u_thrust << endl;
	cout << "u_torque:	 " << u_torque << endl;
	cout << "torqGlob " << (R.mul(u_torque)) << endl;

	
	// Checking that the system equations are the inverse of the endogenous: done
	matrix e3 = {{0},{0},{1}};
	ex accel = (GRAVITY_G * e3 - R * (u_thrust * e3 / mass)).evalm();
	ex angAcc = (J_inertia.inverse() * (u_torque - Omega * J_inertia * omega)).evalm();

	cout << "From inputs u_* back to dynamics: " << endl;
	cout << "accel: " << accel << endl;	
	cout << "angAcc: " << angAcc << endl << endl;


	// checking angular velocity consistency
	cout << "omega in v0: " << state.p << ", " << state.q << ", " << state.r << endl;
	matrix omegaP0 = ex_to<matrix>((R * omega).evalm());
	cout << "omega in p0: " << omegaP0 << endl;
	cout << "omega in v0: " << vectorVrepTransform(omegaP0) << endl;
}
