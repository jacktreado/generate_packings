#include "packing.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

#define DEBUG

#define N 31
#define NDIM 3
#define seed 49

int main(){

	double phi0 = 0.05;
	double dphi = phi0;

	// packing object instantiation
	packing p1(N,NDIM,phi0,seed);

	cout << "starting test jamming md..." << endl;
	// md evoluion
	double t = 400.0;
	double tmp0 = 1.0;
	double dt = 0.1;
	double ep = 10;	
	int plotskip = 500;
	double Utol = N*1e-8;
	double Ktol = N*1e-20;

	// setup simulation
	p1.set_ep(ep);
	p1.set_md_time(dt);
	p1.set_dtmax(10.0);
	p1.set_plotskip(plotskip);
	p1.open_xyz("test.xyz");

	// test for jamming
	p1.jamming_finder(t,dphi,Utol,Ktol);


	return 0;
}