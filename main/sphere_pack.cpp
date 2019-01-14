/*

	Jam N spheres, either bidisperse or monodisperse
	using packing.h 

	BY Jack Treado
	08/27/2018

*/

#include "packing.h"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]){
	cout << "@@ Beginning main for sphere packing..." << endl;
	
	// read in data
	string N_str = argv[1];			// number of particles
	string NDIM_str = argv[2];		// dimension of space
	string phi0_str = argv[3];		// initial packing fraction
	string dphi_str = argv[4];		// delta phi (for growth)
	string alpha_str = argv[5];		// bidispersity ratio (alpha = 1 -> monodisperse)
	string seed_str = argv[6];		// seed
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

	stringstream alphass(alpha_str);
	alphass >> alpha;

	stringstream s1ss(seed_str);
	s1ss >> seed;

	// jamming variables
	int plotskip;
	double ep,dt,Utol,Ktol,t;

	// set parameters
	ep = 10.0;			// energy scale (units of kbt)
	t = 5e3;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	plotskip = 2e3;		// # of steps to skip plotting
	Utol = N*1e-16;		// potential energy tolerance
	Ktol = N*1e-30;		// kinetic energy tolerance

	// NLCL parameters
	int nc,nnu;
	nnu = 100;
	if (N < 1000)
		nc = -1;
	else if (N < 1e5)
		nc = 3;
	else if (N < 1e7)
		nc = 4;
	else
		nc = 5;	

	// packing object instantiation
	cout << "@@ Instantiating packing object..." << endl;
	packing pack(N,NDIM,alpha,phi0,nc,nnu,seed);

	// setup simulation
	pack.set_ep(ep);
	pack.set_md_time(dt);
	pack.set_dtmax(10.0);
	pack.set_plotskip(plotskip);	

	// Bring system to jammed state
	pack.jamming_finder(t,dphi,Utol,Ktol);
	
	// open output files
	pack.open_config(config_str.c_str());
	pack.open_stat(stat_str.c_str());	

	// print output to files
	pack.print_stat();
	pack.print_config();

	cout << "@@ Leaving main for sphere packing!" << endl;
	return 0;
}