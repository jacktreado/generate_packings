// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

int main() {
	// local variables for packing
	string fstr = "../test_input_files/residue_input_test.dat";
	string cfgstr = "residue_cfg_test.dat";
	string statstr = "residue_stat_test.dat";
	string enstr = "residue_test_Energy.dat";
	string xyzstr = "residue_test.xyz";
	int N = 12;
	int dof = 6;
	int nc = -1;
	if (nc >= 80)
		nc = 4;
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
	tmp0 = 0.01;		// initial temperature
	plotskip = 1000;		// # of steps to skip plotting
	phi0 = 0.1;		// initial packing fraction
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
	if (nc > 0) {
		respack.update_neighborlist();
		respack.print_neighborlist();
	}
	respack.rb_jamming_finder(tmp0, NT, dphi, Utol, Ktol);
	// respack.free_fire(tmp0,t);
	// respack.free_md(tmp0,t);

	// output stat and config
	respack.open_stat(statstr.c_str());
	respack.open_config(cfgstr.c_str());

	// print data
	respack.print_stat();
	respack.print_config();

	return 0;
}
