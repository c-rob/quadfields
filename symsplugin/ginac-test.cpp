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
using namespace std;
using namespace GiNaC;

int main()
{
	char fBuffer[128];
	ifstream vectFile("vector-field.txt");
	vectFile.getline(fBuffer, sizeof(fBuffer));
	cout << fBuffer << endl;

	vectFile.close();


	/*
	symbol t("t");
	symtab table;
	table["t"] = t;
	parser reader(table);

	ex e = reader("{{cos(t), -sin(t)}, {sin(t), cos(t)}}");

	ex eNum = e.subs(t == 6);
    cout << "R(1) = " << eNum << endl;	
	double d = ex_to<numeric>(evalf(eNum)).to_double();
    cout << "R(1) = " << d << endl;	

    ex e1 = 2*x*x-4*x+3;
    cout << "e1(7) = " << e1.subs(x == 6) << endl;	

	matrix m = {{x , y}, {pow(x, 2), e1}};
	cout << "m(x,y) = " << m << endl;
	
	exmap varsMap;
	varsMap[x] = 3;
	varsMap[y] = 4;
	cout << "m(3,4) = " << m.subs(varsMap) << endl;
	*/
}
