/************************************************************************
* File parser for a vector field. Each line contains a component of the *
* field in a direction, in the format:                                  *
* V_x                                                                   *
* V_y                                                                   *
* Expression can contain 'x', 'y' variables                             *
************************************************************************/

#include <iostream>
#include <fstream>
#include <ginac/ginac.h>

using namespace GiNaC;
using namespace std;


// Symbolic funcions to evaluate at each interface call
matrix flatOut_D1;
matrix flatOut_D2;
matrix flatOut_D3;
matrix flatOut_D4;


int main()
{

	string line;
	ifstream vectFile("vector-field.txt");
	vector<string> vectFieldStr;
	int nVars;

	// Prepare the GiNaC parser
	symbol x("x");
	symbol y("y");
	symtab table;
	table["x"] = x;
	table["y"] = y;
	parser reader(table);

	// Get the first lines in the file as a vector
	while (getline(vectFile, line)) {
		vectFieldStr.push_back(line);
	}
	vectFile.close();
	nVars = vectFieldStr.size();

	// Fill a symbolic matrix
	matrix vectFieldSym(nVars, 1);
	for (int i = 0; i < nVars; ++i) {
		ex e = reader(vectFieldStr[i]);
		vectFieldSym.set(i, 0, e);
	}


	// Save the first flat output derivative sigma=V(x)
	flatOut_D1 = vectFieldSym;
	 
	// TODO: write a getJacobian() function given the symbolic matrix in input
	
	// debug
	//for (auto l: vectFieldStr) {
	//	cout << l << endl;
	//}
	cout << flatOut_D1 << endl;
	//exmap m;
    //m[x] = -2;
    //m[y] = 6;
	//cout << vectFieldSym.subs(m) << endl;

}
