// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]) {
	// input variables
	int N;
	double Utol;

	// input args
	string N_str = argv[1];			// number of particles
	string input_str = argv[2];		// input file
	string Utol_str = argv[3];		// Utol string

	// get int for number of particles
	stringstream Nss(N_str);
	Nss >> N;

	stringstream Utolss(Utol_str);
	Utolss >> Utol;

	string cfgstr = "dimer_cfg.test";		// config file
	string statstr = "dimer_stat.test";		// stat file
	string enstr = "dimer_en.test";
	string xyzstr = "dimer.xyz";

	// local variables for packing
	int dof = 5;
	int nc = -1;
	int seed = 1;

	// initialize rigid body packing
	rigidbody dimer_pack(input_str,N,dof,nc,seed);

	// set MD parameters
	double ep, dt, tmp0, phi0, dphi, Ktol;
	int plotskip, NT;

	ep = 10.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	tmp0 = 0.001;		// initial temperature
	plotskip = 2e3;		// # of steps to skip plotting
	phi0 = 0.1;		// initial packing fraction
	dphi = 0.001;		// initial packing fraction step
	Ktol = N * 1e-25;
	dimer_pack.rb_scale(phi0);

	// setup simulation
	dimer_pack.set_ep(ep);
	dimer_pack.set_md_time(dt);
	dimer_pack.set_dtmax(10.0);
	dimer_pack.set_plotskip(plotskip);

	// update xW based on xM
	dimer_pack.pos_frot();

	// open xyz and en files
	dimer_pack.open_en(enstr.c_str());
	dimer_pack.open_xyz(xyzstr.c_str());

	// run md
	dimer_pack.rb_jamming_finder(tmp0, NT, dphi, Utol, Ktol);
	// respack.free_fire(tmp0,t);
	// respack.free_md(tmp0,t);

	// output stat and config
	dimer_pack.open_stat(statstr.c_str());
	dimer_pack.open_config(cfgstr.c_str());

	// print data
	dimer_pack.print_stat();
	dimer_pack.print_config();

	return 0;
}
