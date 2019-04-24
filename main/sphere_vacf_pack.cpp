/*

	Jam N spheres, either bidisperse or monodisperse
	using packing.h,
	AND THEN calculate dynamical matrix

	BY Jack Treado
	02/24/2019

*/

#include "packing.h"

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
	string vacfNT_str = argv[7];	// vacf NT run limit
	string vacfT0_str = argv[8];	// vacf T0 run
	string seed_str = argv[9];		// seed
	string config_str = argv[10];	// file to save stats
	string stat_str = argv[11];		// file to save cluster list
	string vacf_str = argv[12];		// file to save dynamical matrix	

	// get numerical values for input variables
	int N,NDIM,seed,vacfNT;
	double dphi,phi0,alpha,dphiDM,vacfT0,vacfNT_tmp;

	// get string stream objects
	stringstream Nss(N_str);
	stringstream NDIMss(NDIM_str);
	stringstream phi0ss(phi0_str);
	stringstream dphiss(dphi_str);
	stringstream alphass(alpha_str);
	stringstream dphiDMss(dphiDM_str);
	stringstream vacfNTss(vacfNT_str);
	stringstream vacfT0ss(vacfT0_str);
	stringstream s1ss(seed_str);

	// put string stream values into variables
	Nss >> N;	
	NDIMss >> NDIM;
	phi0ss >> phi0;
	dphiss >> dphi;
	alphass >> alpha;
	dphiDMss >> dphiDM;
	vacfNTss >> vacfNT_tmp;
	vacfT0ss >> vacfT0;
	s1ss >> seed;

	// cast vacfNT and nsamp to ints
	vacfNT = (int)vacfNT_tmp;

	cout << "vacfNT_tmp = " << vacfNT_tmp << endl;
	cout << "vacfNT = " << vacfNT << endl;

	// jamming variables
	int plotskip,NT,vsave;
	double ep,dt,Utol,Ktol,T0;

	// set parameters
	ep = 1.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.01;			// time step (units of md time)
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

	// open output object
	ofstream dmobj(vacf_str.c_str());
	if (!dmobj.is_open()){
		cout << "file string " << vacf_str << " is not a valid file name, ending" << endl;
		return 1;
	}

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
		pack.dynamical_matrix(dmobj,h);
		cout << "calc complete! printed to " << vacf_str << ". " << endl;

		// // output also the velocity autocorrelation function
		vsave = 64;
		cout << "Calculating VACF for comparison..." << endl;
		cout << "vacfNT = " << vacfNT << endl;
		cout << "vsave = " << vsave << endl;
		pack.calc_vacf(vacfNT,vsave,1,vacfT0,dmobj);
		cout << "VACF calculation complete!" << endl;
	}	


	cout << "@@ Leaving main for sphere packing!" << endl;
	return 0;
}