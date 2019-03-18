// rotate rigidbody test code

#include <iostream>
#include "rigidbody.h"

using namespace std;

int main(){

	// system variables
	int N = 8;
	int DOF = 6;
	string input_str = "/Users/JackTreado/_pv/cluster/rigidbody/io/res_input_N8_seed1.dat";
	string xyzstr = "residue_rotate_test.xyz";

	// load system with object
	rigidbody testpack(input_str,N,DOF,-1,1);

	// open xyz file
	testpack.open_xyz(xyzstr.c_str());

	// output initial W frame positions
	testpack.rigidbody_xyz();

	// update xW based on xM
	testpack.pos_frot();

	// output final W frame positions, check if same
	testpack.rigidbody_xyz();

	// // rotate particle 0 around the horn
	double dtheta = 0.1;
	
	// rotate forward
	testpack.rotate_single_particle_xyz(0,0,dtheta);
	testpack.rotate_single_particle_xyz(1,1,dtheta);
	testpack.rotate_single_particle_xyz(2,2,dtheta);
	testpack.rigidbody_xyz();

	// rotate forward
	testpack.rotate_single_particle_xyz(0,0,dtheta);
	testpack.rotate_single_particle_xyz(1,1,dtheta);
	testpack.rotate_single_particle_xyz(2,2,dtheta);
	testpack.rigidbody_xyz();

	// rotate back
	testpack.rotate_single_particle_xyz(0,0,-dtheta);
	testpack.rotate_single_particle_xyz(1,1,-dtheta);
	testpack.rotate_single_particle_xyz(2,2,-dtheta);
	testpack.rigidbody_xyz();

	// rotate back
	testpack.rotate_single_particle_xyz(0,0,-dtheta);
	testpack.rotate_single_particle_xyz(1,1,-dtheta);
	testpack.rotate_single_particle_xyz(2,2,-dtheta);
	testpack.rigidbody_xyz();


	return 0;
}