// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]) {
	// N
	int N;

	// input args
	string N_str = argv[1];			// number of particles
	string input_str = argv[2];		// input file

	// get int for number of particles
	stringstream Nss(N_str);
	Nss >> N;

	// local variables for packing
	string xyzstr = "residue_test.xyz";
	int dof = 6;
	int nc = -1;
	if (N >= 64)
		nc = 3;
	int seed = 1;

	// initialize rigid body packing
	rigidbody respack(input_str,N,dof,nc,seed);

	// setup energy & viz output
	respack.open_xyz(xyzstr.c_str());

	// update xW based on xM
	respack.pos_frot();

	// output initial W frame positions
	cout << "Printing out initial conditions..." << endl;
	respack.rigidbody_xyz();

	// scramble based on initial seed
	cout << "Scrambling..." << endl;
	respack.scramble();

	// print xyz
	cout << "Printing out final conditions." << endl;
	respack.rigidbody_xyz();

	return 0;
}
