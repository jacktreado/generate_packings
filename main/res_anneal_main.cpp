/*

	Jam N residues, with NLCL option,
	but with an annealing protocol

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
	string fskip_str = argv[5];		// number of annealing steps
	string phimin_str = argv[6];	// minimum packing fraction for annealing
	string tmp0_str = argv[7];		// annealing temperature
	string seed_str = argv[8];		// seed
	string input_str = argv[9];		// file with initial atomic coordinates
	string config_str = argv[10];	// file to save stats
	string stat_str = argv[11];		// file to save cluster list

	// get numerical values for input variables
	int N,NDIM,seed,fskip;
	double dphi,phi0,alpha,phimin,tmp0;

	stringstream Nss(N_str);
	Nss >> N;

	stringstream NDIMss(NDIM_str);
	NDIMss >> NDIM;

	stringstream phi0ss(phi0_str);
	phi0ss >> phi0;

	stringstream dphiss(dphi_str);
	dphiss >> dphi;

	stringstream fskipss(fskip_str);
	fskipss >> fskip;

	stringstream tmp0ss(tmp0_str);
	tmp0ss >> tmp0;

	stringstream phiminss(phimin_str);
	phiminss >> phimin;

	stringstream s1ss(seed_str);
	s1ss >> seed;

	// output input parameters
	cout << endl;
	cout << "@@ INPUT PARAMETERS:" << endl;
	cout << "N = " << N << endl;
	cout << "NDIM = " << NDIM << endl;
	cout << "phi0 = " << phi0 << endl;
	cout << "dphi = " << dphi << endl;
	cout << "fskip = " << fskip << endl;
	cout << "tmp0 = " << tmp0 << endl;
	cout << "phimin = " << phimin << endl;
	cout << "seed = " << seed << endl;
	cout << "@@" << endl << endl;

	// jamming variables
	int plotskip,nc,dof,NT;
	double ep,dt,Utol,Ktol;	

	// set parameters
	dof = 6;			// number of degrees of freedom per particle
	ep = 10.0;			// energy scale (units of kbt)
	NT = 1e8;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	plotskip = 1e3;		// # of steps to skip plotting
	Utol = N*1e-16;		// potential energy tolerance
	Ktol = N*1e-30;		// kinetic energy tolerance

	// rigid body object instantiation
	cout << "@@ Instantiating packing object..." << endl;

	// NLCL parameters, if system is large enough
	if (N < 300)
		nc = -1;	
	else if (N < 1000)
		nc = 3;
	else if (N < 5000)
		nc = 4;
	else
		nc = 5;
	
	// instantiate rigidbody md object
	rigidbody rp(input_str,N,dof,nc,seed);

	// setup simulation	
	rp.rb_scale(phi0);
	rp.set_ep(ep);
	rp.set_md_time(dt);
	rp.set_dtmax(10.0);
	rp.set_plotskip(plotskip);
	rp.pos_frot();		

	// Bring system to jammed state
	rp.rb_anneal(tmp0,NT,fskip,phimin,dphi,Utol,Ktol);
	
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
