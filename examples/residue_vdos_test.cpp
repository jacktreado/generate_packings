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

	// dynamical matrix variables
	string dm_str = "/Users/JackTreado/_pv/vdos/mixed/res_dm.dat"; 		// dynamical matrix file
	double dphiDM = 1e-14;

	// get int for number of particles
	stringstream Nss(N_str);
	Nss >> N;

	// local variables for packing
	string cfgstr = "residue_cfg.test";
	string statstr = "residue_stat.test";
	string enstr = "residue_Energy.test";
	string xyzstr = "residue_test.xyz";
	int dof = 6;
	int nc = -1;
	if (N >= 64)
		nc = 3;
	int seed = 1;

	// initialize rigid body packing
	rigidbody respack(input_str,N,dof,nc,seed);

	// Do short md to check if code is set up correctly

	// set MD parameters
	double ep, dt, tmp0, phi0, dphi, Utol, Ktol;
	int plotskip, NT;

	ep = 10.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.025;			// time step (units of md time)
	tmp0 = 0.01;		// initial temperature
	plotskip = 1000;	// # of steps to skip plotting
	phi0 = 0.1;			// initial packing fraction
	dphi = 0.001;		// initial packing fraction step
	Utol = N * 1e-20;
	Ktol = N * 1e-30;	

	// scale particles
	respack.rb_scale(phi0);

	// setup simulation
	respack.set_ep(ep);
	respack.set_md_time(dt);	
	respack.set_dtmax(10.0); // ADD FROM res_pack_nlcl
	respack.set_plotskip(plotskip);		

	// setup energy & viz output
	respack.open_en(enstr.c_str());
	respack.open_xyz(xyzstr.c_str());

	// output initial W frame positions
	respack.rigidbody_xyz();

	// update xW based on xM
	respack.pos_frot();

	// output final W frame positions, check if same
	respack.rigidbody_xyz();

	// run md
	respack.rb_jamming_finder(tmp0, NT, dphi, Utol, Ktol);

	// print config, dynamical matrix
	int isjammed = 0;
	isjammed = respack.get_isjammed();
	double h = 1e-8;


	if (isjammed == 1){
			// open output files
		respack.open_config(cfgstr.c_str());
		respack.open_stat(statstr.c_str());	

		// print output to files
		respack.print_stat();
		respack.print_config();

		// scale system slightly
		phi0 = respack.get_phi();
		cout << "scaling phi from phi = " << phi0 << " to phi = " << phi0 + dphiDM << endl;
		respack.rb_scale(phi0 + dphiDM);

		int NTmax = 5e5;
		cout << "minimizing potential energy..." << endl;
		respack.rb_fire_umin(NTmax,Ktol);
		cout << "energy minimization complete!" << endl;

		cout << "starting to calculate the dynamical matrix..." << endl;
		respack.rb_dynamical_matrix(dm_str,h);
		cout << "calc complete! printed to " << dm_str << ". " << endl;	
	}



	cout << "@@ Leaving main for residue packing!" << endl;
	return 0;
}