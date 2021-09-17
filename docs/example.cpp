#include <stdio.h>
#include <tnt/tnt.h>
#include <flt/system.hpp>
#include <flt/derivatives.hpp>

using namespace std;
using namespace TNT;
using namespace FLT;

/** \example example.cpp
 *
 * Obtains the output of a closed loop fuzzy system and its jacobian matrix.
 *
 * <i>Compile</i>:
 *    g++ example.cpp ../src/membership.cpp ../src/rule.cpp ../src/system.cpp ../src/utilities.cpp ../src/derivatives.cpp -o example
 * 
 * <i>Use</i>:
 *    example plant.txt controller.txt x1 x2 ... xn ('n' is the order of the plant)
 * 
 * <i>Test</i>:
 *    ./example ../matlab_utilities/Plant.txt ../matlab_utilities/Controller.txt 0 0
 * 
 * <i>Output</i>:\n
 *    <i>Closed loop output</i>:\n
 *    0.0409887\n
 *    0.340054
 * 
 *    <i>Closed loop Jacobian matrix</i>:\n
 *    -2.58935  1.98356\n
 *    2.30224  -7.94545
*/

int main(int argc, char **argv)
{	
	if (argc < 4)
	{
		cout << "Error. Some input aguments are missing.\n\n";
		cout << "example plant.txt controller.txt x1 x2 ... xn\n";
		cout << "\t where 'n' is the order of the plant (number of outputs)" << endl;
		return 1;
	}
	
	char *plant = argv[1];
	System P = TXT2System(plant);
	size_t n = P.outputs();
	if (!n)
	{
		cout << "Error reading the plant file." << endl;
		return 1;
	}
	
	if (argc != 3+n)
	{
		cout << "Error. Some input aguments are missing.\n\n";
		cout << "example plant.txt controller.txt x1 x2 ... xn\n";
		cout << "\t where 'n' is the order of the plant (number of outputs)" << endl;
		return 1;
	}
	
	char *controller = argv[2];
	System C = TXT2System(controller);
	size_t m = C.outputs();
	if (!m)
	{
		cout << "Error reading the controller file." << endl;
		return 1;
	}
	
	if (P.inputs() != n+m)
	{
		cout << E_NoCoherent << endl;
		return 1;
	}
	
	Array1D<double> X(n);
	for (size_t i=0; i<n; i++)
		X[i] = atof(argv[i+3]);
	
	Array1D<double> dX = evaluate(&X[0], P, C);
	
	cout << "\nClosed loop output:\n";
	for (size_t i=0; i<n; i++)
		cout << dX[i] << "\n";
	
	Array2D<double> J = jacobian(P, C, &X[0]);
	cout << "\nClosed loop Jacobian matrix:\n";
	for (size_t i=0; i<n; i++)
	{
		for (size_t j=0; j<n; j++)
			cout << J[i][j] << "  ";
		cout << "\n";
	}
	cout << endl;
	
	return 0;
}
