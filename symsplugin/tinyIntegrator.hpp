/****************************************************************************
* This small class wraps in one place the code for the basic euler          *
* integration. That is an useful operation for debugging. The internal type *
* is GiNaC's 'numeric' type. Use it:                                        *
*                                                                           *
*     TinyIntegrator variable("name", 0.005); // 5 ms                       *
*     variable.setInitialState(matrix{{2},{3},{0}});                        *
*     for (;;)                                                              *
*         variable.update(velocityVect);                                    *
*     cout << variable << endl;                                             *
*                                                                           *
****************************************************************************/

#include <iostream>
#include <string>
#include <ginac/ginac.h>

using std::string;
using GiNaC::ex;

class TinyIntegrator {

	private:
		string name;
		ex val;
		double dt;
		bool initialized;

		TinyIntegrator(string varName, double timeStep, ex initialState, bool initDone):
			name(varName), val(initialState), dt(timeStep), initialized(initDone) {}

		friend std::ostream& operator<<(std::ostream& os, const TinyIntegrator& var);
		
	public:

		TinyIntegrator(string varName, double timeStep):
			TinyIntegrator(varName, timeStep, 0, false) {}

		TinyIntegrator(string varName, double timeStep, ex initialState):
			TinyIntegrator(varName, timeStep, initialState, false) {}

		ex get(void) const {
			return val;
		}

		GiNaC::matrix getMat(void) const {
			return GiNaC::ex_to<GiNaC::matrix>(val);
		}

		string getName(void) const {
			return name;
		}

		void setInitialState(ex initalState) {
			val = initalState;
			initialized = true;
		}

		ex update(ex rate) {
			val = (val + rate * dt).evalm();
			return val;
		}

		bool isInitialized() {
			return initialized;
		}
};
