/*

	Jam monodisperse spheres 
	using packing.h 

	BY Jack Treado
	07/02/2018

*/

#include "packing.h"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#define N 500
#define NDIM 3
#define seed 1
#define NCL 3

int main(int argc, char *argv[]){
	cout << "@@ Beginning nl printing test for monodisperse packing..." << endl;
	
	string config_str = "config_nl.dat";	// file to save stats
	string stat_str = "stat_nl.dat";		// file to save cluster list
	string xyz_str = "nl.xyz";

	// main variables
	double phi0 = 0.01;
	double dphi = 0.005;
	double alpha = 1.4;

	// jamming variables
	int plotskip;
	double ep,dt,Utol,Ktol,t;

	// set MD parameters
	ep = 10.0;			// energy scale (units of kbt)
	t = 5.0;			// total amount of time (units of sim time)
	dt = 0.05;			// time step (units of md time)
	plotskip = 200;		// # of steps to skip plotting
	Utol = N*1e-8;		// potential energy tolerance
	Ktol = N*1e-20;		// kinetic energy tolerance

	// neighbor list variables
	int nc = 4;
	int nnu = 10;
	double rcut = 3.5;

	// packing object instantiation
	packing p1(N,NDIM,alpha,phi0,nc,nnu,rcut,seed);

	// setup simulation
	p1.set_ep(ep);
	p1.set_md_time(dt);
	p1.set_dtmax(10.0);
	p1.set_plotskip(plotskip);	
	p1.open_xyz(xyz_str.c_str());

	cout << "Printing starting info" << endl;
	p1.print_xyz();

	cout << "getting cell neighbors" << endl;
	p1.get_cell_neighbors();

	cout << "getting cell positions" << endl;
	p1.get_cell_positions();
	cout << "printing cell positions" << endl;
	p1.print_cell_pos();

	cout << "update cell" << endl;
	p1.update_cell();
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
	p1.jamming_finder_nn(t,dphi,Utol,Ktol);

	// open output files
	p1.open_config(config_str.c_str());
	p1.open_stat(stat_str.c_str());	

	// print output to files
	p1.print_stat();
	p1.print_config();

	cout << "@@ Leaving nl printing test for monodisperse packing!" << endl;
	return 0;
}