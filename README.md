# SymDeriv V-REP plugin

This is a V-REP plugin for small symbolic computations written in c++
It implements equations 8-11 of paper:
Vector Field Following for Quadrotors using Differential Flatness
Dingjiang Zhou and Mac Schwager

## Dependencies:
	* libginac
	* libcln

## Linux installation:
	* Install the dependencies
	* Set in the Makefile the variable VREPDIR to the V-REP installation
	    directory
	* cd to the symsplugin directory
	* Run "make mklib" to build the library or "make install" to build
		and copy to the v-rep installation directory

### Notes:
	* The plugin assumes that the quadrotor is a single dynamic shape
        in the vrep scene that is called "Quadrirotor".
	* See "Field controller" scene for usage
	* The plugin has been written for V-REP 3.4.0. A newer API might require
		minor changes in the code.
	* Euler angles singularity at 90° pitch
