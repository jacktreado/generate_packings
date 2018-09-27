// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

int main(){
	// local variables for packing
	string fstr = "/Users/JackTreado/_pv/cluster/res/io/res_input_N128_seed1.dat";
	int N = 128;
	int dof = 6;
	int nc = 4;
	int seed = 1;

	// initialize rigid body packing
	rigidbody respack(fstr,N,dof,nc,seed);

	// output info to see if done right	
	cout << "Printing all variables from packing class: " << endl;
	respack.print_vars();

	// cout << "\nPrinting rotational variables: " << endl;
	// respack.print_rot();

	cout << "\nPrinting cell info: " << endl;
	respack.print_cell();

	cout << "\n Printing cell pos: " << endl;
	respack.print_cell_pos();

	cout << "\nPrinting clabel: " << endl;
	respack.print_clabel();

	return 0;
}