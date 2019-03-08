/*

	Script to calculate dynamical matrix from an existing 
	jammed packing of spheres

	1. Read in file
	2. Scale packing fraction slightly
	3. Minimize potential energy
	4. calculate dynamical matrix
	5. IF LOOP: rescale phi, minimize potential energy, calculate again

*/

#include <iostream>
#include <iomanip>
#include "packing.h"
#include <string>
#include <sstream>

using namespace std;

#define NDIM 3

int main(int argc, char* argv[]){
	// load in command line arguments
	string Nstr = argv[1];			// number of particles
	string dphistr = argv[2];		// delta phi value
	string cfgstr = argv[3];		// config string
	string dmstr = argv[4];			// dynamical matrix output file

	// open test energy file too
	string enstr = "sphere_energy.test";	

	// parse string
	stringstream Nss(Nstr);
	stringstream dphiss(dphistr);

	// initialize main parameters
	int N, seed, plotskip, NTmax;
	double phi0, dphi, h, Ktol, ep, dt;

	Nss >> N;
	dphiss >> dphi;

	// set parameters
	ep = 10.0;			// energy scale (units of kbt)
	dt = 0.025;			// time step (units of md time)
	Ktol = N*1e-30;		// kinetic energy tolerance
	seed = 1;			// initial seed
	h = 1e-10;			// small change parameter
	plotskip = 1e2;		// plotskip	
	NTmax = 5e5; 		// number of energy minimize steps to try

	// instantiate object
	cout << "instantiating object" << endl;
	packing jammed_sphere_packing(cfgstr,NDIM,seed);
	cout << "object instantiation complete" << endl;

	// setup simulation
	jammed_sphere_packing.set_ep(ep);
	jammed_sphere_packing.set_md_time(dt);
	jammed_sphere_packing.set_dtmax(10.0);
	jammed_sphere_packing.set_plotskip(plotskip);	

	// open energy obj
	jammed_sphere_packing.open_en(enstr);

	// scale packing fraction	
	phi0 = jammed_sphere_packing.get_phi();
	cout << "scaling phi from phi = " << phi0 << " to phi = " << phi0 + dphi << endl;
	jammed_sphere_packing.scale_sys(dphi);

	// energy minimizer
	// cout << "set alpha0" << endl;
	// jammed_sphere_packing.set_alpha0(0.75);
	// jammed_sphere_packing.set_fdec(0.9999);
	cout << "minimizing potential energy..." << endl;
	jammed_sphere_packing.fire_umin(NTmax,Ktol);
	cout << "energy minimization complete!" << endl;

	// vibrational density of states
	cout << "starting to calcualte the dynamical matrix..." << endl;
	jammed_sphere_packing.dynamical_matrix(dmstr,h);
	cout << "calc complete! printed to " << dmstr << ", now ending program." << endl;	

	return 0;
}
