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
	vector <symbol> vars = {x, y};
	nVars = vars.size();
	table["x"] = x;
	table["y"] = y;
	parser reader(table);

	// Get the first lines in the file as a vector
	while (getline(vectFile, line)) {
		vectFieldStr.push_back(line);
	}
	vectFile.close();

	// Fill a symbolic matrix
	matrix vectFieldSym(nVars, 1);
	for (int i = 0; i < nVars; ++i) {
		ex e = reader(vectFieldStr[i]);
		vectFieldSym.set(i, 0, e);
	}

	// Save the first flat output derivative d(sigma)/dt=V(x)
	flatOut_D1 = vectFieldSym;

	// Compute next derivatives
	genNextDerivative(vars, flatOut_D1, flatOut_D2);
	genNextDerivative(vars, flatOut_D2, flatOut_D3);
	genNextDerivative(vars, flatOut_D3, flatOut_D4);

	
	cout << flatOut_D1 << endl;
	cout << flatOut_D2 << endl;
	cout << flatOut_D3 << endl;
	cout << flatOut_D4 << endl;
}

