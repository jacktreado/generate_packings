// test for reading in info to rigid body class

#include "rigidbody.h"
#include "vec3.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]) {
	// input variables
	int N;
	double h,dphiDM;

	// input args
	string N_str = argv[1];			// number of particles
	string input_str = argv[2];		// input file
	string h_str = argv[3];
	string dphi_str = argv[4];

	// dynamical matrix variables
	string dm_str = "/Users/JackTreado/_pv/vdos/mixed/dimer_dm.dat";

	// get int for number of particles
	stringstream Nss(N_str);
	Nss >> N;

	stringstream hss(h_str);
	hss >> h;

	stringstream dphiss(dphi_str);
	dphiss >> dphiDM;

	// local variables for packing
	string cfgstr = "dimer_cfg.test";
	string statstr = "dimer_stat.test";
	string enstr = "dimer_Energy.test";
	string xyzstr = "dimer_test.xyz";
	int dof = 5;
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
	phi0 = 0.4;			// initial packing fraction
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

	// run md
	respack.rb_jamming_finder(tmp0, NT, dphi, Utol, Ktol);

	// print config, dynamical matrix
	int isjammed = 0;
	isjammed = respack.get_isjammed();


	if (isjammed == 1){
			// open output files
		respack.open_config(cfgstr.c_str());
		respack.open_stat(statstr.c_str());	

		// print output to files
		respack.print_stat();
		respack.print_config();

		// scale system slightly
		phi0 = respack.get_phi();
		cout << "scaling phi from phi = " << phi0 << " to phi = " << phi0 + dphiDM << " with dphiDM = " << dphiDM << endl;
		respack.rb_scale(phi0 + dphiDM);

		int NTmax = 5e5;
		cout << "minimizing potential energy..." << endl;
		respack.rb_fire_umin(NTmax,Ktol);
		cout << "energy minimization complete!" << endl;

		cout << "starting to calculate the dynamical matrix..." << endl;
		respack.open_xyz(xyzstr.c_str());
		respack.dimer_dynamical_matrix(dm_str,h);
		cout << "calc complete! printed to " << dm_str << ". " << endl;	
	}



	cout << "@@ Leaving main for dimer packing VDOS calc!" << endl;
	return 0;
}