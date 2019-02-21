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
#include "packing"
#include <string>
#include <>
#include <sstream>

using namespace std;

#define NDIM 3

int main(int argc, char* argv[]){
	// load in command line arguments
	string Nstr = argv[1];				// number of particles
	string dphistr = argv[2];			// delta phi value
	string cfgstr = argv[3];			// config string
	string outputptrn = argv[4];		// pattern for output 

	// parse string
	stringstream Nss(Nstr);
	stringstream dphiss(dphistr);

	// initialize main parameters
	int N;
	double dphi;

	Nss >> N;
	dphiss >> dphi;

	// instantiate object
	packing jammed_pack(cfgstr,N,NDIM,dphi);


}
