// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]) {
	// main variables
	int N = 1500;
	string input_str = "/Users/JackTreado/Jamming/ProteinVoids/cluster/rigidbody/io/dimer_input_N1500_seed1.dat";
	string cfgstr = "dimer_nlcl_cfg.test";		// config file
	string statstr = "dimer_nlcl_stat.test";		// stat file
	string enstr = "dimer_nlcl_en.test";
	string xyzstr = "dimer_nlcl.xyz";

	// local variables for packing
	int dof = 5;
	int nc = 4;
	int seed = 1;

	// initialize rigid body packing
	rigidbody dimer_pack(input_str,N,dof,nc,seed);

	// set MD parameters
	double ep, dt, tmp0, phi0, dphi, Ktol, Utol;
	int plotskip, NT;

	ep = 1.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.01;			// time step (units of md time)
	tmp0 = 1e-2;		// initial temperature
	plotskip = 1e2;		// # of steps to skip plotting
	phi0 = 0.62;		// initial packing fraction
	dphi = 0.001;		// initial packing fraction step
	Utol = N * 1e-8;
	Ktol = N * 1e-30;

	cout << "IN MAIN: updating phi from " << dimer_pack.get_phi() << " to phi0 = " << phi0 << endl;
	dimer_pack.rb_scale(phi0);
	cout << "IN MAIN: new phi = " << dimer_pack.get_phi() << endl;

	// setup simulation
	dimer_pack.set_ep(ep);
	dimer_pack.set_md_time(dt);
	dimer_pack.set_dtmax(10.0);
	dimer_pack.set_plotskip(plotskip);

	// update xW based on xM
	dimer_pack.pos_frot();	

	// run md
	// dimer_pack.rb_jamming_finder(tmp0, NT, dphi, Utol, Ktol);
	cout << "IN MAIN: Running fire to relax U" << endl;
	dimer_pack.free_fire(tmp0, Utol, 1e2);

	cout << "IN MAIN: Running MD to check E con" << endl;

	// open xyz and en files
	dimer_pack.open_en(enstr.c_str());
	dimer_pack.open_xyz(xyzstr.c_str());
	dimer_pack.free_md(tmp0, 5e4);

	// output stat and config
	dimer_pack.open_stat(statstr.c_str());
	dimer_pack.open_config(cfgstr.c_str());

	// print data
	dimer_pack.print_stat();
	dimer_pack.print_config();

	return 0;
}
