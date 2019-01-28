// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]){
	string ntstr = argv[1];
	string t0str = argv[2];

	// print strings
	cout << "NT input = " << ntstr << endl;
	cout << "T0 input = " << t0str << endl;

	// read in using string stream
	stringstream ntss(ntstr);
	stringstream t0ss(t0str);	

	// read in inputs
	double NT_tmp;
	double T0;
	ntss >> NT_tmp;
	t0ss >> T0;

	cout << "NT val = " << NT_tmp << endl;
	cout << "T0 val = " << T0 << endl;

	// cast double NT to int
	int NT;
	NT = (int)NT_tmp;
	cout <<"int NT val = " << NT << endl;

	return 0;

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
	double ep,dt,phi0,dphi,Utol,Ktol;
	int plotskip,nnu;

	ep = 10.0;			// energy scale (units of kbt)
	dt = 0.05;			// time step (units of md time)
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
