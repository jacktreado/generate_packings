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
	string cfg_str = argv[3];		// config file
	string stat_str = argv[4];		// stat file
	string Utol_str = argv[5];		// Utol string

	// get int for number of particles
	stringstream Nss(N_str);
	Nss >> N;

	stringstream Utolss(Utol_str);
	Utolss >> Utol;

	// local variables for packing
	int dof = 5;
	int nc = -1;
	int seed = 1;

	// initialize rigid body packing
	rigidbody respack(input_str,N,dof,nc,seed);

	// set MD parameters
	double ep, dt, tmp0, phi0, dphi, Ktol;
	int plotskip, NT;

	ep = 10.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	tmp0 = 0.001;		// initial temperature
	plotskip = 2e4;		// # of steps to skip plotting
	phi0 = 0.1;		// initial packing fraction
	dphi = 0.005;		// initial packing fraction step
	Ktol = N * 1e-20;
	respack.rb_scale(phi0);

	// setup simulation
	respack.set_ep(ep);
	respack.set_md_time(dt);
	respack.set_dtmax(10.0);
	respack.set_plotskip(plotskip);

	// update xW based on xM
	respack.pos_frot();

	// run md
	respack.rb_jamming_finder(tmp0, NT, dphi, Utol, Ktol);
	// respack.free_fire(tmp0,t);
	// respack.free_md(tmp0,t);

	// output stat and config
	respack.open_stat(stat_str.c_str());
	respack.open_config(cfg_str.c_str());

	// print data
	respack.print_stat();
	respack.print_config();

	return 0;
}
