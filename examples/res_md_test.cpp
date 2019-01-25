// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

int main(){
	// local variables for packing
	string fstr = "/Users/JackTreado/_pv/cluster/rigidbody/rcp/config/res_rcp_config_N16_seed1.dat";
	string enstr = "res_md_en.test";
	string xyzstr = "res_md_test.xyz";
	int N = 16;
	int dof = 6;
	int nc = -1;
	if (nc >= 64)
		nc = 3;
	int seed = 1;

	// initialize rigid body packing
	rigidbody respack(fstr,N,dof,nc,seed);

	// Do short md to check if code is set up correctly

	// set MD parameters
	double ep,dt,tmp0,phi0,dphi,Utol,Ktol;
	int plotskip,NT,nnu;

	ep = 10.0;			// energy scale (units of kbt)
	NT = 2e5;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	tmp0 = 1e-4;		// initial temperature
	plotskip = 500;		// # of steps to skip plotting
	nnu = 5;

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

	return 0;
}
