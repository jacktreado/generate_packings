/*

	Jam N spheres, either bidisperse or monodisperse
	using packing.h 

	BY Jack Treado
	08/27/2018

*/


#include "rigidbody.h"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]){
	cout << "@@ Beginning main for res packing..." << endl;
	
	// read in data
	string N_str = argv[1];			// number of particles
	string NDIM_str = argv[2];		// dimension of space
	string phi0_str = argv[3];		// initial packing fraction
	string dphi_str = argv[4];		// delta phi (for growth)
	string seed_str = argv[5];		// seed
	string input_str = argv[6];		// file with initial atomic coordinates
	string config_str = argv[7];	// file to save stats
	string stat_str = argv[8];		// file to save cluster list

	// get numerical values for input variables
	int N,NDIM,seed;
	double dphi,phi0,alpha;

	stringstream Nss(N_str);
	Nss >> N;

	stringstream NDIMss(NDIM_str);
	NDIMss >> NDIM;

	stringstream phi0ss(phi0_str);
	phi0ss >> phi0;

	stringstream dphiss(dphi_str);
	dphiss >> dphi;

	stringstream s1ss(seed_str);
	s1ss >> seed;

	// jamming variables
	int plotskip,nc,dof,NT;
	double ep,dt,Utol,Ktol,tmp0;	

	// set parameters
	dof = 6;			// number of degrees of freedom per particle
	ep = 10.0;			// energy scale (units of kbt)
	NT = 1e8;			// total amount of time (units of sim time)
	tmp0 = 0.01;		// initial temperature
	dt = 0.05;			// time step (units of md time)
	plotskip = 1e3;		// # of steps to skip plotting
	Utol = N*1e-8;		// potential energy tolerance
	Ktol = N*1e-20;		// kinetic energy tolerance

	// rigid body object instantiation
	cout << "@@ Instantiating packing object..." << endl;

	// NLCL parameters, if system is large enough
	if (N >= 80){
		cout << "@@ System large enough to warrant NLCL..." << endl;
		if (N < 200)
			nc = 4;			// number of cells along one direction
		else
			nc = 5;

	}
	else
		nc = -1;
	
	// instantiate rigidbody md object
	rigidbody rp(input_str,N,dof,nc,seed);

	// setup simulation	
	rp.rb_scale(phi0);
	rp.set_ep(ep);
	rp.set_md_time(dt);
	rp.set_dtmax(10.0);
	rp.set_plotskip(plotskip);		

	// Bring system to jammed state
	if (nc > 0){
		rp.update_neighborlist();
		rp.print_neighborlist();
	}
	rp.rb_jamming_finder(tmp0,NT,dphi,Utol,Ktol);
	
	// open output files
	int isjammed = 0;
	isjammed = rp.get_isjammed();
	if (isjammed == 1){
		rp.open_config(config_str.c_str());
		rp.open_stat(stat_str.c_str());	

		// print output to files
		rp.print_stat();
		rp.print_config();
	}
	else{
		cout << "isjammed = 0, packing did not find jammed state!" << endl;

	cout << "@@ Leaving main for res packing!" << endl;
	return 0;
	}
}