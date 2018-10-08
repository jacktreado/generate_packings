// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

int main() {
	// local variables for packing
  	string fstr = "/Users/JackTreado/Jamming/ProteinVoids/cluster/res/io/res_input_N48_seed1.dat";
	string cfgstr = "residue_cfg.test";
	string statstr = "residue_stat.test";
	string enstr = "residue_Energy.test";
	string xyzstr = "residue_test.xyz";
	int N = 48;
	int dof = 6;
	int nc = 4;
	if (N >= 40)
		nc = 3;
	int seed = 1;

	// initialize rigid body packing
	rigidbody respack(fstr, N, dof, nc, seed);

	// Do short md to check if code is set up correctly

	// set MD parameters
	double ep, dt, tmp0, phi0, dphi, Utol, Ktol;
	int plotskip, NT;

	ep = 10.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	tmp0 = 0.001;		// initial temperature
	plotskip = 100;		// # of steps to skip plotting
	phi0 = 0.5;		// initial packing fraction
	dphi = 0.001;		// initial packing fraction step
	Utol = N * 1e-8;
	Ktol = N * 1e-20;
	respack.rb_scale(phi0);

	// setup simulation
	respack.set_ep(ep);
	respack.set_md_time(dt);
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
