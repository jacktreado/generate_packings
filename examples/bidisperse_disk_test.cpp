/*

	Jam N bidisperse disks in 2 dimensions

	BY Jack Treado
	05/16/2019

*/

#include "packing.h"

using namespace std;

int main(){
	cout << "beginning of main for bidisperse disk packing." << endl;

	// get numerical values for input variables
	int N,NDIM,seed,nc,nnu;
	double dphi,phi0,alpha;
	string config_str, xyz_str;

	// parameters
	N 		= 32;
	NDIM 	= 2;
	phi0 	= 0.01;
	dphi 	= 1e-3;
	alpha 	= 1.4;
	seed 	= 14234123;

	// strings
	config_str = "disk_cfg.test";
	xyz_str = "disk.xyz";

	// jamming variables
	int plotskip,NT;
	double ep,dt,Utol,Ktol,T0;

	// set parameters
	ep = 1.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.01;			// time step (units of md time)
	plotskip = 1e3;		// # of steps to skip plotting
	Utol = N*1e-16;		// potential energy tolerance
	Ktol = N*1e-30;		// kinetic energy tolerance

	// NLCL parameters
	nc = -1;
	nnu = 1;

	// packing object instantiation
	cout << "Instantiating packing object..." << endl;
	packing pack(N,NDIM,alpha,phi0,nc,nnu,seed);

	// setup simulation
	pack.set_ep(ep);
	pack.set_md_time(dt);
	pack.set_dtmax(10.0);
	pack.set_plotskip(plotskip);
	pack.open_xyz(xyz_str.c_str());

	// Bring system to jammed state
	pack.jamming_finder(NT,dphi,Utol,Ktol);

	// check if jammed
	int isjammed = 0;
	isjammed = pack.get_isjammed();

	// if jammed print
	if (isjammed){
		cout << "opening configuration file " << config_str << " to print stats and configuration" << endl;
		pack.open_config(config_str.c_str());

		cout << "printing data to configuration file" << endl;
		pack.print_data();
	}
	else
		cout << "Jamming configuration not found for input seed = " << seed << endl;

	// end main
	cout << "end of main for bidisperse disk packing." << endl;
	return 0;
}