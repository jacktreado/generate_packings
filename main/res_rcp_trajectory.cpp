/*

	file to do res rigid body MD
	and output trajectory file in .xyz format

*/

#include "rigidbody.h"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;


int main(int argc, char *argv[]){
	cout << "@@ BEGINNING main for res md trajectory ... " << endl;

	// local variables
	int N,NT,seed,dof,nc,nnu,plotskip;
	double NT_tmp,T0,ep,dt;

	// read in options
	string N_str = argv[1];				// NUMBER OF PARTICLES
	string NT_str = argv[2];			// MAX # OF TIME STEPS
	string T0_str = argv[3];			// TEMPERATURE
	string seed_str = argv[4];			// VELOCITY SEED
	string input_str = argv[5];			// INPUT CFG FILE
	string xyz_str = argv[6];			// TRAJECTORY FILE
	string final_cfg_str = argv[7];		// FINAL CONFIG FILE (EQUILIBRATED)

	// parse numeric options using stringstream
	stringstream Nss(N_str);
	stringstream NTss(NT_str);
	stringstream T0ss(T0_str);
	stringstream seedss(seed_str);

	// stream in values
	Nss >> N;	
	NTss >> NT_tmp;	
	T0ss >> T0;	
	seedss >> seed;

	// cast NT value to integer, in case of scientific notation input
	NT = (int)NT_tmp;

	// set parameters to initialize rigidbody packing
	dof = 6;
	nc = -1;
	if (N >= 64)
		nc = 3;

	// instantiate rigidbody packing object
	rigidbody trajobj(input_str,N,dof,nc,seed);

	// open xyz file
	trajobj.open_xyz(xyz_str.c_str());

	// set parameters to initialize MD
	ep = 10.0;			// energy scale (units of kbt)
	dt = 0.05;			// time step (units of md time)
	plotskip = 500;		// # of steps to skip plotting
	nnu = 5;			// NLCL update if needed

	// expected MB USAGE
	cout << "@@ Expect " << 0.002*(NT/plotskip)*N << " MB for the trajectory file..." << endl;

	// setup simulation
	trajobj.set_ep(ep);
	trajobj.set_md_time(dt);
	trajobj.set_plotskip(plotskip);	

	// output initial W frame positions
	trajobj.rigidbody_xyz();

	// update xW based on xM
	trajobj.pos_frot();

	// output final W frame positions, check if same
	trajobj.rigidbody_xyz();

	// run free md	
	if (nc > 0){
		trajobj.update_neighborlist();
		trajobj.print_neighborlist();
	}
	trajobj.free_md(T0,NT,nnu);

	// print final config to config file
	trajobj.open_config(final_cfg_str.c_str());
	trajobj.print_config();

	cout << "@@ ENDING main for res md trajectory ... " << endl;
	return 0;
}
