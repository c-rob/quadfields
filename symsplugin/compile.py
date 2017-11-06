#!/usr/bin/env python3
"""
Run the compiler for a vrep-plugin with a given name:
	compile.sh <name.cpp>
"""

import sys
import os

if __name__ == '__main__':
	
	fileName = sys.argv[1]
	exeName = os.path.splitext(fileName)[0]

	os.system('gcc ' + fileName + ' -o ' + exeName + ' -lcln -lstdc++ -lginac')
