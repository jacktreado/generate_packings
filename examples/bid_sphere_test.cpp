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

int main(){
	cout << "@@ Beginning main for sphere packing..." << endl;
	
	// save data
	string config_str = "sphere_cfg.test";	// file to save stats
	string stat_str = "sphere_stat.test";		// file to save cluster list
	string xyz_str = "sphere.xyz";
	string en_str = "sphere_energy.test";

	// get numerical values for input variables
	int N,NDIM,seed;
	double dphi,phi0,alpha;

	N = 15;
	NDIM = 3;
	seed = 1;
	dphi = 0.001;
	phi0 = 0.05;
	alpha = 1.26;

	// jamming variables
	int plotskip;
	double ep,dt,Utol,Ktol,t;

	// set parameters
	ep = 10.0;			// energy scale (units of kbt)
	t = 5e3;			// total amount of time (units of sim time)
	dt = 0.01;			// time step (units of md time)
	plotskip = 500;		// # of steps to skip plotting
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
	pack.scale_sys(dphi);

	// open xyz and energy files
	pack.open_xyz(xyz_str);
	pack.open_en(en_str);

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