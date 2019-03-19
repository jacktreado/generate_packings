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
	double ep, dt, tmp0, phi0, dphi, Utol, Ktol, phimin;
	int plotskip, NT, fskip;

	ep = 1.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.01;			// time step (units of md time)	
	plotskip = 500;	// # of steps to skip plotting
	phi0 = 0.01;			// initial packing fraction
	dphi = 0.001;		// initial packing fraction step
	Utol = N * 1e-16;
	Ktol = N * 1e-30;	

	// ANNEALING PARAMETERS
	tmp0 = 1e-16;		// initial temperature (and 10*kick temperature)
	fskip = 20;			// number of steps between fire minimizations
	phimin = 0.3;		// minimum packing fraction to try to anneal

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
	respack.rb_anneal(tmp0, NT, fskip, phimin, dphi, Utol, Ktol);
	// respack.free_fire(tmp0,t);
	// respack.free_md(tmp0,t);

	// output stat and config
	respack.open_stat(statstr.c_str());
	respack.open_config(cfgstr.c_str());

	// print data
	respack.print_stat();
	respack.print_config();

	// print xyz
	respack.rigidbody_xyz();

	return 0;
}
