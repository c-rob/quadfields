# SymDeriv V-REP plugin

This is a V-REP plugin for small symbolic computations written in c++
It implements equations 8-11 of paper:
Vector Field Following for Quadrotors using Differential Flatness
Dingjiang Zhou and Mac Schwager

## Dependencies:
* libginac-dev
* libcln-dev

## Linux installation:
* Install the dependencies
* Set in the Makefile the variable VREPDIR to your V-REP installation
	directory
* cd to the symsplugin directory
* Run "make mklib" to build the library or "make install" to build
	and copy to the v-rep installation directory

### Notes:
* The plugin assumes that the quadrotor is a single dynamic shape
        in the vrep scene that is called "Quadrirotor".
* See "Field controller" scene for usage
* To use the scene, also modify the path of the file vector-field.txt inside
	the script in the V-REP scene
* The plugin has been written for V-REP 3.4.0. A newer API might require
		minor changes in the code.
* Euler angles singularity at 90° pitch
