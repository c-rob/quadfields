
#include "tinyIntegrator.hpp"

using std::cout;
using std::endl;

std::ostream& operator<<(std::ostream& os, const TinyIntegrator& var)  
{  
	os << var.getName() << ": " << var.get();  
	return os;  
}


/*
int main() {


	TinyIntegrator x("pippo", 0.01, GiNaC::matrix({{1},{2}}));

	cout << x << endl;
	x.update(GiNaC::matrix({{1.3},{-4}}));
	x.update(GiNaC::matrix({{1.3},{-4}}));
	cout << x << endl;
}
*/
