/*

	Jam N spheres, either bidisperse or monodisperse
	using packing.h,
	AND THEN calculate dynamical matrix

	BY Jack Treado
	02/24/2019

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
	string dphiDM_str = argv[6];	// increases system size by dphi to calc DM
	string seed_str = argv[7];		// seed
	string config_str = argv[8];	// file to save stats
	string stat_str = argv[9];		// file to save cluster list
	string dm_str = argv[10];		// file to save dynamical matrix

	// get numerical values for input variables
	int N,NDIM,seed;
	double dphi,phi0,alpha,dphiDM;

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

	stringstream dphiDMss(dphiDM_str);
	dphiDMss >> dphiDM;

	stringstream s1ss(seed_str);
	s1ss >> seed;

	// jamming variables
	int plotskip,NT,vsave;
	double ep,dt,Utol,Ktol,T0;

	// set parameters
	ep = 10.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.025;			// time step (units of md time)
	plotskip = 5e3;		// # of steps to skip plotting
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
	pack.jamming_finder(NT,dphi,Utol,Ktol);

	// print config, dynamical matrix
	int isjammed = 0;
	isjammed = pack.get_isjammed();
	double h = 1e-8;	

	// print to analytical form as well, just to check
	string dm2_str = "/Users/JackTreado/Jamming/ProteinVoids/vdos/mixed/analytical_check.dat";
	string vacf_str = "/Users/JackTreado/Jamming/ProteinVoids/vdos/mixed/sphere_vacf.dat";
	ofstream vacf_obj(vacf_str.c_str());

	if (isjammed == 1){
			// open output files
		pack.open_config(config_str.c_str());
		pack.open_stat(stat_str.c_str());	

		// print output to files
		pack.print_stat();
		pack.print_config();

		// scale system slightly
		phi0 = pack.get_phi();
		cout << "scaling phi from phi = " << phi0 << " to phi = " << phi0 + dphiDM << endl;
		pack.scale_sys(dphiDM);

		int NTmax = 5e5;
		cout << "minimizing potential energy..." << endl;
		pack.fire_umin(NTmax,Ktol);
		cout << "energy minimization complete!" << endl;

		cout << "starting to calculate the dynamical matrix..." << endl;
		pack.dynamical_matrix(dm_str,h);
		cout << "calc complete! printed to " << dm_str << ". " << endl;

		// analytical DM, to compare
		cout << "calculating DM using analytical form..." << endl;
		pack.analytical_dm(dm2_str);
		cout << "calc complete! printed to " << dm2_str << ", now ending program." << endl;

		// output also the velocity autocorrelation function
		NT = 163840;
		T0 = 1e-20;
		vsave = NT/32768;	// make sure a power of 2 for FFT
		cout << "Calculating VACF for comparison..." << endl;	
		pack.calc_vacf(NT,vsave,T0,vacf_obj);
		cout << "VACF calculation complete!" << endl;
	}	


	cout << "@@ Leaving main for sphere packing!" << endl;
	return 0;
}