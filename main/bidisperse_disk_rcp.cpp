/*

	Jam N bidisperse disks in 2 dimensions

	USING ATHERMAL COMPRESSION

	BY Jack Treado
	05/16/2019

*/

#include "packing.h"

// define NDIM to be 2, always dealing with bidisperse disks
#define NDIM 2

using namespace std;

int main(int argc, char *argv[]){
	// read in data
	string N_str = argv[1];			// number of particles
	string phi0_str = argv[2];		// initial packing fraction
	string dphi_str = argv[3];		// delta phi (for growth)
	string alpha_str = argv[4];		// radii ratio (should be != 1)
	string seed_str = argv[5];		// seed
	string config_str = argv[6];	// file to save stats and configuration

	// get numerical values for input variables
	int N,seed;
	double dphi,phi0,alpha;

	stringstream Nss(N_str);
	Nss >> N;

	stringstream phi0ss(phi0_str);
	phi0ss >> phi0;

	stringstream dphiss(dphi_str);
	dphiss >> dphi;

	stringstream alphass(alpha_str);
	alphass >> alpha;

	stringstream s1ss(seed_str);
	s1ss >> seed;

	// print something to console
	cout << endl << endl << endl << endl;
	cout << "===================================" << endl << endl << endl;
	cout << " RUNNING N = " << N << " BIDISPERSE DISKS" << endl;
	cout << " WITH SEED = " << seed << endl << endl << endl;
	cout << "===================================" << endl << endl << endl;
	cout << endl << endl << endl << endl;

	// jamming variables
	int plotskip,NT;
	double ep,dt,Utol,Ktol,T0;

	// set parameters
	ep = 1.0;			// energy scale (units of kbt)
	NT = 5e6;			// total amount of time (units of sim time)
	dt = 0.01;			// time step (units of md time)
	plotskip = 5e3;		// # of steps to skip plotting
	Utol = N*1e-16;		// potential energy tolerance
	Ktol = N*1e-30;		// kinetic energy tolerance

	// packing object instantiation
	cout << "Instantiating packing object..." << endl;
	packing pack(N,NDIM,alpha,phi0,-1,1,seed);

	// setup simulation
	pack.set_ep(ep);
	pack.set_md_time(dt);
	pack.set_dtmax(10.0);
	pack.set_plotskip(plotskip);

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
	cout << "end of main for bidisperse disk packing using athermal compression." << endl;
	return 0;
}