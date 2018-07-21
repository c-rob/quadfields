# FieldFollow V-REP plugin

This is a V-REP plugin that implements a controller for a quadrotor model, written in C++.

The controller is described in:
D. Zhou and M. Schwager, “Vector Field Following
for Quadrotors using Differential Flatness”, In Proc. of the International
Conference on Robotics and Automation (ICRA 14), June, 2014, pp.  6567-6572.

This software is not related in any way with the authors of this paper.

## Dependencies:
* libginac-dev
* libcln-dev

## Linux installation:
* Install the dependencies
* Set in the Makefile the variable VREPDIR to your V-REP installation
	directory
* cd to the symsplugin directory
* Run "make mklib" to build the library or "make install" to build
	and copy it to the V-REP installation directory

### Notes:
* See the scene "field_controller.ttt" for an example of usage.
* To use the scene, also modify the path of the file vector-field.txt inside
	the main quadrotor child script.
* The plugin has been written for V-REP 3.4.0.
* Euler angles singularity at 90° pitch
* Tested with ODE and Bullet <=2.83 (not working with Vortex)
