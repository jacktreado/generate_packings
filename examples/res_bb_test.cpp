#include "backbone.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

int main (int argc, char *argv[]){
	// opening print statement
	cout << "starting bb test main!" << endl;

	// number of particles
	int N = 3;

	// input file path	
	string inputstr = "/Users/JackTreado/_pv/backbone/bb_io/res_bb_input_N3_seed1.dat";

	// config and stat strings
	string cfgstr = "bb_config.test";
	string statstr = "bb_stat.test";
	string xyzstr = "bb.xyz";
	string enstr = "bb_en.test";

	// simulation parameters
	int dof = 6;
	int nc = -1;
	int seed = 1;

	// instantiate packing object
	backbone bb_pack(inputstr,N,dof,nc,seed);

	// open xyz string
	bb_pack.open_xyz(xyzstr);
	bb_pack.open_en(enstr);

	// set MD parameters
	double ep, dt, tmp0, phi0, dphi, Ktol, kbl, kba, kda, falpha0, falpha1;
	int plotskip, NT;

	ep = 10.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	tmp0 = 0.001;		// initial temperature
	plotskip = 10;		// # of steps to skip plotting
	phi0 = 0.01;			// initial packing fraction
	dphi = 0.005;		// initial packing fraction step
	Ktol = N * 1e-20;
	falpha0 = 0.99;
	falpha1 = 0.99;

	cout << "set spring constants" << endl;
	kbl = 1;
	kba = 1;
	kda = 1;	
	bb_pack.set_kbl(kbl);
	bb_pack.set_kba(kba);
	bb_pack.set_kda(kda);

	cout << "set initial falpha to be " << falpha0 << endl;
	bb_pack.set_falpha(falpha0);

	cout << "scaling system to initial packing fraction phi0 = " << phi0 << endl;
	bb_pack.rb_scale(phi0);

	// setup simulation
	cout << "setting MD parameters" << endl;
	bb_pack.set_ep(ep);
	bb_pack.set_md_time(dt);
	bb_pack.set_dtmax(10.0);
	bb_pack.set_plotskip(plotskip);

	// update xW based on xM
	cout << "rotating system to xM for MD purposes" << endl;
	bb_pack.pos_frot();

	// run md
	cout << "relaxing topology of backbone" << endl;
	bb_pack.top_relax();

	// open stat and config files
	cout << "opening stat and config files" << endl;
	bb_pack.open_stat(statstr.c_str());
	bb_pack.open_config(cfgstr.c_str());


	// print data
	cout << "printing stat and config files" << endl;
	bb_pack.print_stat();
	bb_pack.print_config();

	cout << "ending bb test main!" << endl;
	return 0;
}