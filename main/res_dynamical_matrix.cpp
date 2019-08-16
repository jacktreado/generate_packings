/*

	Main file to calculate dynamical matrix of a given residue packing.

	Uses numerical approximation for the Hessian of the potential energy, 
	from https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm

	Output eigenvalues and eigenvectors to [FILE NAME].dm
	
*/

// include files
#include "rigidbody.h"
#include "vec3.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

// namespace
using namespace std;

// simulation parameters
const int dof 		= 6;		// degree of freedom per particle
const int NT 		= 1e6;		// number of time steps
const int plotskip 	= 5e3;		// number of steps to skip during simulation
const int nc 		= -1;		// never use NLCL list (for now)

const double ep 	= 1.0;		// potential energy scale
const double dt 	= 0.025;	// time step
const double tmp0 	= 1e-4;		// temperature (does not effect result)
const double phi0 	= 0.4;		// initial packing fraction
const double dphi 	= 0.001;	// packing fraction increment

int main(int argc, char *argv[]) {
	// input variables
	int N, seed;
	double h, dphiDM, phiJ;

	// input args
	string N_str 		= argv[1];			// number of particles
	string h_str 		= argv[2];			// pertubation scale
	string dphi_str 	= argv[3];			// DM phi increment scale
	string seed_str 	= argv[4];			// initial seed for packing
	string input_str 	= argv[5];			// string for input file
	string cfg_str 		= argv[6];			// string for final configuration
	string dm_str 		= argv[7];			// string for dynamical matrix file

	// get int for number of particles
	stringstream Nss(N_str);
	stringstream hss(h_str);
	stringstream dphiss(dphi_str);
	stringstream seedss(seed_str);

	hss >> h;
	Nss >> N;
	dphiss >> dphiDM;
	seedss >> seed;

	// initialize rigid body packing
	rigidbody respack(input_str,N,dof,nc,seed);

	// energy tolerances
	double Utol, Ktol;
	Utol = N * 1e-20;
	Ktol = N * 1e-30;	

	// scale particles
	respack.rb_scale(phi0);

	// setup simulation
	respack.set_ep(ep);
	respack.set_md_time(dt);	
	respack.set_dtmax(10.0); // ADD FROM res_pack_nlcl
	respack.set_plotskip(plotskip);		

	// update xW based on xM
	respack.pos_frot();

	// run md to find initial jammed state
	respack.rb_jamming_finder(tmp0, NT, dphi, Utol, Ktol);

	// print config, dynamical matrix
	int isjammed = 0;
	isjammed = respack.get_isjammed();

	if (isjammed == 1){
		// open output files
		respack.open_config(cfg_str.c_str());

		// print output to files
		respack.print_config();

		// scale system slightly
		phiJ = respack.get_phi();
		cout << "scaling phi from phiJ = " << phiJ << " to phi = " << phiJ + dphiDM << " with dphiDM = " << dphiDM << endl;
		respack.rb_scale(phiJ + dphiDM);

		int NTmax = 5e5;
		cout << "minimizing potential energy..." << endl;
		respack.rb_fire_umin(NTmax,Ktol);
		cout << "energy minimization complete!" << endl;

		cout << "starting to calculate the dynamical matrix..." << endl;
		respack.rb_dynamical_matrix(dm_str,h);
		cout << "calc complete! printed to " << dm_str << ". " << endl;	
	}

	cout << "@@ Leaving main for residue packing!" << endl;
	return 0;
}