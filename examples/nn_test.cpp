/*

	Jam spheres (mono- or bi-disperse)
	using packing.h 

	BY Jack Treado
	07/02/2018

	updated:
	09/27/18

*/

#include "packing.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

#define N 100
#define NDIM 3
#define seed 3

int main(int argc, char *argv[]){
	cout << "@@ Beginning nl printing test for sphere packing..." << endl;
	
	string config_str = "config_nl.test";	// file to save stats
	string stat_str = "stat_nl.test";		// file to save cluster list
	string xyz_str = "nl.xyz";

	// main variables
	double phi0 = 0.01;
	double dphi = 0.001;
	double alpha = pow(2,1.0/3.0);

	// jamming variables
	int plotskip;
	double ep,dt,Utol,Ktol,t;

	// set MD parameters
	ep = 10.0;			// energy scale (units of kbt)
	t = 5.0;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	plotskip = 1000;		// # of steps to skip plotting
	Utol = N*1e-8;		// potential energy tolerance
	Ktol = N*1e-20;		// kinetic energy tolerance

	// neighbor list variables
	int nc = 3;
	int nnu = 100;

	// packing object instantiation
	cout << "Instiating object" << endl;
	cout << "alpha = " << alpha << endl;
	packing p1(N,NDIM,alpha,phi0,nc,nnu,seed);

	// setup simulation
	cout << "Setting ep" << endl;
	p1.set_ep(ep);

	cout << "Setting md time" << endl;
	p1.set_md_time(dt);

	cout << "Setting dtmax for FIRE" << endl;
	p1.set_dtmax(10.0);

	cout << "Setting plotskip variable" << endl;
	p1.set_plotskip(plotskip);	

	cout << "Opening XYZ File" << endl;
	p1.open_xyz(xyz_str.c_str());

	cout << "Printing starting info to xyz file" << endl;
	p1.print_xyz();

	cout << "printing cell positions" << endl;
	p1.print_cell_pos();

	cout << "printing" << endl;
	p1.print_cell();
	p1.print_clabel();
	p1.print_celln();	
	p1.print_cell_neighbors();	

	cout << "update neighbor list" << endl;
	p1.update_neighborlist();

	cout << "printing" << endl;	
	p1.print_neighborlist();	
	p1.print_nl_xyz();

	// Bring system to jammed state
	cout << "Beginnning jamming" << endl;
	p1.jamming_finder(t,dphi,Utol,Ktol);

	// open output files
	int isjammed = 0;
	isjammed = p1.get_isjammed();
	if (isjammed == 1){
		p1.open_config(config_str.c_str());
		p1.open_stat(stat_str.c_str());	

		// print output to files
		p1.print_stat();
		p1.print_config();
	}
	else
		cout << "packing did not find jammed state! " << endl;

	cout << "@@ Leaving nl printing test for sphere packing!" << endl;
	return 0;
}