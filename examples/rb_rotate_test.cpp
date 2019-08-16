// rotate rigidbody test code

#include <iostream>
#include "rigidbody.h"

using namespace std;

int main(){

	// system variables
	int N = 8;
	int DOF = 6;
	string input_str = "/Users/JackTreado/_pv/cluster/rigidbody/io/dimer_input_N3_seed1.dat";
	string xyzstr = "dimer_rotate_test.xyz";

	// load system with object
	rigidbody testpack(input_str,N,DOF,-1,1);

	// open xyz file
	testpack.open_xyz(xyzstr.c_str());

	// update xW based on xM
	testpack.pos_frot();

	// // rotate particle 0 around the horn
	double dtheta = 0.3;

	// plot original
	testpack.rigidbody_xyz();
	testpack.rotate_single_particle(0,1,dtheta);
	testpack.rigidbody_xyz();
	testpack.rotate_single_particle(0,1,-dtheta);
	testpack.rigidbody_xyz();

	// rotate multiple times
	for (int i=0; i<40; i++){
		testpack.rotate_single_particle(0,1,dtheta);
		testpack.rigidbody_xyz();
	}
	for (int i=0; i<40; i++){
		testpack.rotate_single_particle(1,2,dtheta);
		testpack.rigidbody_xyz();
	}
	for (int i=0; i<40; i++){
		testpack.rotate_single_particle(2,3,dtheta);
		testpack.rigidbody_xyz();
	}


	return 0;
}