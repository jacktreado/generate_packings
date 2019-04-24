// test for reading in info to rigid body class

#include "rigidbody.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

#define NDIM 3
#define DOF 6

using namespace std;

int main(int argc, char *argv[]) {
	// input args
	string N_str = argv[1];			// number of particles
	string phi0_str = argv[2];		// initial packing fraction
	string dphi_str = argv[3];		// delta phi (for growth)
	string dphiDM_str = argv[4];	// increases system size by dphi to calc DM
	string vacfNT_str = argv[5];	// vacf NT run limit
	string vacfT0_str = argv[6];	// vacf T0 run
	string tsave_str = argv[7];		// MD steps between save
	string seed_str = argv[8];		// seed
	string input_str = argv[9];		// input string
	string config_str = argv[10];	// file to save jammed config
	string stat_str = argv[11];		// file to save states about final config
	string vel_str = argv[12];		// file to save velocity information

	// get numerical values for input variables
	int N,seed,vacfNT,tsave;
	double dphi,phi0,dphiDM,vacfT0,vacfNT_tmp;

	// get string stream objects
	stringstream Nss(N_str);
	stringstream phi0ss(phi0_str);
	stringstream dphiss(dphi_str);
	stringstream dphiDMss(dphiDM_str);
	stringstream vacfNTss(vacfNT_str);
	stringstream vacfT0ss(vacfT0_str);
	stringstream tsavess(tsave_str);
	stringstream s1ss(seed_str);

	// put string stream values into variables
	Nss >> N;
	phi0ss >> phi0;
	dphiss >> dphi;
	dphiDMss >> dphiDM;
	vacfNTss >> vacfNT_tmp;
	vacfT0ss >> vacfT0;
	tsavess >> tsave;
	s1ss >> seed;

	// cast vacfNT and nsamp to ints
	vacfNT = (int)vacfNT_tmp;

	// initialize rigid body packing (-1 turns neighbor list off )
	rigidbody respack(input_str,N,DOF,-1,seed);

	// Do short md to check if code is set up correctly

	// set MD parameters
	double ep, dt, tmp0, Utol, Ktol;
	int plotskip, NT;

	ep = 1.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of packing time (MD steps)
	dt = 0.01;			// time step (units of md time)
	tmp0 = 0.01;		// initial temperature
	plotskip = 5000;	// # of steps to skip plotting
	phi0 = 0.1;			// initial packing fraction
	dphi = 0.001;		// initial packing fraction step
	Utol = N * 1e-20;
	Ktol = N * 1e-30;	

	// scale particles
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

	// print config, dynamical matrix
	int isjammed = 0;
	isjammed = respack.get_isjammed();
	double h = 1e-8;


	if (isjammed == 1){
		// open output files
		respack.open_config(config_str.c_str());
		respack.open_stat(stat_str.c_str());	

		// print output to files
		respack.print_stat();
		respack.print_config();

		// open velocity output file
		ofstream mdobj(vel_str.c_str());
		if (mdobj.fail()){
			cout << "MD object failed to open with string " << vel_str << endl;
			exit(1);
		}

		// scale system slightly
		phi0 = respack.get_phi();
		cout << "scaling phi from phi = " << phi0 << " to phi = " << phi0 + dphiDM << endl;
		respack.rb_scale(phi0 + dphiDM);

		int NTmax = 5e5;
		cout << "minimizing potential energy..." << endl;
		respack.rb_fire_umin(NTmax,Ktol);
		cout << "energy minimization complete!" << endl;

		
		// run MD trajectory for VACF calculation
		cout << "running MD at increased packing fraction for NT = " << vacfNT << " steps" << endl;	
		respack.rb_md_velocity(vacfT0,vacfNT,tsave,mdobj);
		cout << "finished running MD on residues for VACF" << endl;


		// close md output file
		mdobj.close();
	}



	cout << "@@ Leaving main" << endl;
	return 0;
}