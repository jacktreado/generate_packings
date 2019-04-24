/*

	Jam N spheres, either bidisperse or monodisperse
	using packing.h,
	AND THEN calculate dynamical matrix

	AND perturb single mode, see if you can recreate VACF

	BY Jack Treado
	04/15/2019

*/

#include "packing.h"

using namespace std;

int main(){
	cout << "@@ Beginning example main to calc sphere vacf from single mode..." << endl;

	// get numerical values for input variables
	int N,NDIM,seed,vacfNT,mode,nc,nnu,vsave;
	double dphi,phi0,alpha,dphiDM,vacfT0;
	string config_str,stat_str,vacf_str;

	// parameters
	N 		= 20;
	NDIM 	= 3;
	phi0 	= 0.1;
	dphi 	= 1e-3;
	alpha 	= 1.26;
	dphiDM  = 1e-4;
	vacfNT 	= (int)pow(2,18);
	vacfT0 	= 1e-20;
	seed 	= 2;
	vsave 	= 64;

	// nontrivial mode to perturb
	mode 	= 3;

	// strings
	config_str = "sphere_cfg.test";
	stat_str = "sphere_stat.test";
	if (mode > -1)
		vacf_str = "/Users/JackTreado/_pv/vdos/mixed/single_mode_pert.dat";
	else
		vacf_str = "/Users/JackTreado/_pv/vdos/mixed/sphere_vacf_test.dat";

	// jamming variables
	int plotskip,NT;
	double ep,dt,Utol,Ktol,T0;

	// set parameters
	ep = 1.0;			// energy scale (units of kbt)
	NT = 5e8;			// total amount of time (units of sim time)
	dt = 0.01;			// time step (units of md time)
	plotskip = 5e3;		// # of steps to skip plotting
	Utol = N*1e-16;		// potential energy tolerance
	Ktol = N*1e-30;		// kinetic energy tolerance

	// NLCL parameters
	nc = -1;
	nnu = 1;

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

		pack.set_plotskip(10*plotskip);
		if (mode > -1){
			cout << "perturbing a single mode: mode = " << mode << endl;			
			pack.single_mode_perturbation(dmobj,mode,vacfNT,vsave,vacfT0);
			cout << "Finishing single mode perturbation " << endl;
		}
		else{
			cout << "perturbing all modes...." << endl;			
			pack.all_mode_perturbation(dmobj,vacfNT,vsave,vacfT0);
		}
		
	}	


	cout << "@@ Leaving example main for sphere packing!" << endl;
	return 0;
}