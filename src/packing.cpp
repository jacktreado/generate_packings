/*

	Methods implementation 
	for packing class

	BY Jack Treado

*/

#include "packing.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

#define DEBUG_OFF

const double PI = 3.1415926;


/* 
==================================

	CONSTRUCTORS & DESTRUCTORS		 

================================== 
*/

// initialize packing for rigid body stuff
packing::packing(int n, int dof, int nc, int s){
	// set known member variables
	N = n;
	seed = s;
	NDIM = 3;
	NCN = 26;
	DOF = dof;
	NCL = nc;
	NCELLS = pow(nc,NDIM);
	nnupdate = 1;	
	this->initialize_NC();

	// local variables
	int i,j,d;

	// set other member variables to be changed
	phi = -1;
	rcut = -1;

	// initialize box to -1
	this->initialize_box(-1.0);

	// initialize particles
	this->initialize_particles();

	// setup cell if appropriate
	if (nc < 0)
		this->freeze_nlcl();
	else
		this->setup_nlcl();

	// set std sim variables (energy, plotting, etc.)
	this->setup_std_sim();

	// set std FIRE variables
	this->setup_std_FIRE();

	// set default plot value to 1
	plotit = 1;
}

// N particles, from file (need NDIM)
packing::packing(string &str, int ndim, int s){
	// local variables
	int i,j,d;
	double rad;

	// set parameters to member variables
	NDIM = ndim;
	DOF = NDIM;
	seed = s;
	this->initialize_NC();

	// read in particle positions
	this->read_spheres(str);

	// initial dt = 1
	dt = 1.0;

	// random initial positions & velocities, set force to 0, r to phi0 val
	v = new double*[N];
	F = new double*[N];
	aold = new double*[N];
	pc = new int[N];	
	for (i=0; i<N; i++){
		v[i] = new double[NDIM];
		F[i] = new double[NDIM];
		aold[i] = new double[NDIM];
		pc[i] = 0;
		for (d=0; d<NDIM; d++){
			v[i][d] = drand48();
			F[i][d] = 0.0;
			aold[i][d] = 0.0;
		}
	}

	// initialize energy
	ep = 1.0;
	U = 0.0;
	K = 0.0;

	// initialize contact matrix
	c = new double[NC];
	for (i=0; i<NC; i++)
		c[i] = 0;

	// initialize number of rattlers (not checked at t = 0)
	nr = 0;

	// FIRE variables
	alpha0 = 0.1;
	alpha = alpha0;
	finc = 1.1;
	fdec = 0.5;
	falpha = 0.99;
	np = 0;

	// output variables
	plotskip = 1000;
	isjammed = 0;

	// CL/NL variables: set to nullptr, except for NL, make everything a neighbor to everything
	rcut = -1;
	NCL = -1;
	NCELLS = -1;
	NCN = -1;
	nnupdate = -1;
	g = nullptr;
	clabel = nullptr;
	cellpos = nullptr;
	cellneighbors = nullptr;
	celln = nullptr;
	cell = nullptr;

	neighborlist = new vector<int>[N];
	for (i=0; i<N; i++){
		for (j=i+1; j<N; j++)
			neighborlist[i].push_back(j);
	}

	// set default plot value to 1
	plotit = 1;
}


// N particles using neighbor list
packing::packing(int n, int ndim, double alpha, double phi0, int nc, int nnu, int seed){
	// set initial seed
	srand48(seed);

	// local variables
	int i,j,d;
	double rad,msum,vol,dphi,rc;

	// set parameters to member variables
	N = n;
	NDIM = ndim;
	DOF = NDIM;
	NCL = nc;
	NCELLS = pow(nc,NDIM);
	nnupdate = nnu;
	this->initialize_NC();

	// get rad, # of neighbors, rcut magnitude (set to 2.5x large particle radius)
	if (NDIM == 2){
		rad = sqrt(phi0/(N*PI));
		if (NCL == -1)
			NCN = -1;
		else{
			rc = 2.5;
			NCN = 8;
		}
	}
	else if (NDIM == 3){
		rad = pow((3.0*phi0)/(N*4.0*PI),1.0/NDIM);
		if (NCL == -1)
			NCN = -1;
		else{
			rc = 2.5;
			NCN = 26;
		}
	}		
	else{
		cout << "NDIM = " << NDIM << " not yet supported...." << endl;
		cout << "ENDING PROGRAM.\n";
		throw "NDIM not supported\n";
	}
	cout << "rad = " << rad << endl;

	// get rcut
	rcut = rc*alpha*rad;

	// initial dt = 1
	dt = 1.0;

	// initial box = unit length
	this->initialize_box(1.0);

	// random initial positions & velocities, set force to 0, r to phi0 val
	this->initialize_particles(seed,rad,alpha);	

	// either "freeze" cell list by setting ptrs to null, or allocate memory
	if (NCL == -1)
		this->freeze_nlcl();
	else{
		this->setup_nlcl();
		this->get_cell_neighbors();
		this->get_cell_positions();
		this->update_cell();
		this->update_neighborlist();
	}

	// set std sim variables (energy, plotting, etc.)
	this->setup_std_sim();

	// set std FIRE variables
	this->setup_std_FIRE();

	// set default plot value to 1
	plotit = 1;

	// scale particles to desired size
	vol = 1.0;
	msum = 0;
	for(d=0; d<NDIM; d++)
		vol *= L[d];
	for(i=0; i<N; i++)
		msum += m[i];
	phi = msum/vol;
	dphi = phi0-phi;
	this->scale_sys(dphi);
}






packing::~packing(){
	// delete 1d arrays
	delete [] L;
	delete [] r;
	delete [] m;
	delete [] c;
	delete [] pc;

	// delete 2d arrays
	for (int i=0; i<N; i++){
		delete [] x[i];
		delete [] v[i];
		delete [] F[i];
		neighborlist[i].clear();		
	}
	delete [] x;
	delete [] v;
	delete [] F;

	// delete cell list/neighbor list variables
	delete [] g;
	delete [] celln;
	delete [] clabel;
	for (int i=0; i<NCELLS; i++){
		delete [] cellpos[i];
		delete [] cellneighbors[i];
		cell[i].clear();		
	}
	delete [] cellpos;
	delete [] cellneighbors;
	delete [] cell;
	delete [] neighborlist;	

	xyzobj.close();
	configobj.close();
	statobj.close();
	enobj.close();
}


/* 
==================================
 
		  INITIALIZATION		 

================================== 
*/


void packing::freeze_nlcl(){
	// local variables
	int i,j;

	// CL/NL variables: set to nullptr, except for NL, make everything a neighbor to everything
	rcut = -1;
	NCL = -1;
	NCELLS = -1;
	NCN = -1;
	nnupdate = -1;
	g = nullptr;
	clabel = nullptr;
	cellpos = nullptr;
	cellneighbors = nullptr;
	celln = nullptr;
	cell = nullptr;

	neighborlist = new vector<int>[N];
	for (i=0; i<N; i++){
		for (j=i+1; j<N; j++)
			neighborlist[i].push_back(j);
	}
}

void packing::setup_nlcl(){
	// local variables
	int i,d;
	
	// initialize cell list
	neighborlist = new vector<int>[N];
	cellpos = new double*[N];
	cellneighbors = new int*[NCELLS];
	cell = new vector<int>[NCELLS];
	celln = new int[NCELLS];
	g = new double[NDIM];
	clabel = new int[N];

	for (i=0; i<N; i++)
		clabel[i] = -1;

	for (i=0; i<NDIM; i++)
		g[i] = L[i]/NCL;

	for (i=0; i<NCELLS; i++){
		celln[i] = 1;
		cellpos[i] = new double[NDIM];
		for (d=0; d<NDIM; d++)
			cellpos[i][d] = -1;

		cellneighbors[i] = new int[NCN];
		for (d=0; d<NCN; d++)
			cellneighbors[i][d] = -1;
	}
}

void packing::initialize_particles(){
	// local variables
	int i,d;

	// random initial positions & velocities, set force to 0, r to phi0 val
	x = new double*[N];
	v = new double*[N];
	F = new double*[N];
	aold = new double*[N];
	r = new double[N];
	m = new double[N];
	pc = new int[N];			
	for (i=0; i<N; i++){
		x[i] = new double[NDIM];
		v[i] = new double[NDIM];
		F[i] = new double[NDIM];
		aold[i] = new double[NDIM];
		r[i] = 0.0;	
		m[i] = 0.0;
		pc[i] = 0;
		for (d=0; d<NDIM; d++){
			x[i][d] = L[d]*drand48();
			v[i][d] = 0.0;
			F[i][d] = 0.0;
			aold[i][d] = 0.0;
		}
	}

	// initialize contact matrix
	c = new double[NC];
	for (i=0; i<NC; i++)
		c[i] = 0;
}

void packing::initialize_particles(int seed, double rad, double alpha){
	// set seed for positions
	srand48(seed);

	// local variables
	int i,d;

	// random initial positions & velocities, set force to 0, r to phi0 val
	x = new double*[N];
	v = new double*[N];
	F = new double*[N];
	aold = new double*[N];
	r = new double[N];
	m = new double[N];
	pc = new int[N];			
	for (i=0; i<N; i++){
		x[i] = new double[NDIM];
		v[i] = new double[NDIM];
		F[i] = new double[NDIM];
		aold[i] = new double[NDIM];		
		pc[i] = 0;
		for (d=0; d<NDIM; d++){
			x[i][d] = L[d]*drand48();
			v[i][d] = 0.0;
			F[i][d] = 0.0;
			aold[i][d] = 0.0;
		}

		if (i < round(0.5*N))
			r[i] = rad;	
		else
			r[i] = alpha*rad;
		m[i] = (4.0/3.0)*PI*pow(r[i],3);

		cout << "r[" << i << "] = " << r[i] << ", m[" << i << "] = " << m[i] << endl;
	}

	// initialize contact matrix
	c = new double[NC];
	for (i=0; i<NC; i++)
		c[i] = 0;

	// check mean mass
	cout << "mean mass in initialize_particles = ";
	double mean_mass = this->get_mean_mass();
	cout << mean_mass << endl;
}

void packing::initialize_box(double val){
	// local variables
	int d;

	// initial box = unit length
	L = new double[NDIM];
	for (d=0; d<NDIM; d++)
		L[d] = val;
}

void packing::setup_std_FIRE(){
	// set std FIRE variables
	alpha0 = 0.1;
	alpha = alpha0;
	finc = 1.1;
	fdec = 0.5;
	falpha = 0.99;
	np = 0;
	dtmax = 10*dt;
}

void packing::setup_std_sim(){
	ep = 1.0;
	U = 0.0;
	K = 0.0;
	nr = 0;	
	plotskip = 1000;
	isjammed = 0;
	dt = 1.0;	
}


/* 
==================================
 
			 SETTERS		 

================================== 
*/

void packing::set_rand_c(double p){
	int i,j,ci,cj;
	double r;

	for(i=0; i<N; i++){
		ci = i*N - ((i+1)*(i+2))/2;
		for (j=i+1; j<N; j++){
			cj = ci + j;		
			r = drand48();
			if (r < p)
				c[cj] = 1;
			else
				c[cj] = 0;
		}
	}
}


void packing::set_md_time(double dt0){
	int i,d;
	double m_avg,r_avg,ks,w;

	// get average m
	m_avg = this->get_mean_mass();
	r_avg = this->get_mean_rad();
	ks = ep/(r_avg*r_avg);

	w = sqrt(ks/m_avg);
	dt = dt0/w;

	cout << "* mean m = " << m_avg << endl;
	cout << "* mean r = " << r_avg << endl;
	cout << "* char. spring = " << ks << endl;
	cout << "* char. freq. = " << w << endl;
	cout << "* MD time step dt = " << dt << endl;
}





/* 
==================================

		  	 GETTERS		 

================================== 
*/

double packing::get_distance(int p1, int p2){
	int d;
	double dr,h;

	h = 0;
	dr = 0;
	for (d=0; d<NDIM; d++){
		dr = x[p2][d]-x[p1][d];
		dr -= L[d]*round(dr/L[d]);
		h += dr*dr;
	}

	h = sqrt(h);
	return h;
}

double packing::get_distance(int p1, int p2, double xij[]){
	int d;
	double dr,h;

	h = 0;
	dr = 0;
	for (d=0; d<NDIM; d++){
		dr = x[p2][d]-x[p1][d];
		dr -= L[d]*round(dr/L[d]);
		xij[d] = dr;
		h += dr*dr;
	}

	h = sqrt(h);
	return h;
}

double packing::get_mean_vel(){
	double val = 0.0;
	int i,d;

	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			val += v[i][d];
		}
	}
	val /= N;
	return val;
}

double packing::get_com(){
	double val = 0.0;
	int i,d;

	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			val += x[i][d];
		}
	}
	val /= N;
	return val;
}

double packing::get_mean_mass(){
	int i;
	double val = 0.0;

	for (i=0; i<N; i++)
		val += m[i];

	val /= N;
	return val;
}

double packing::get_mean_rad(){
	int i;
	double val = 0.0;

	for (i=0; i<N; i++)
		val += r[i];

	val /= N;
	return val;
}

double packing::get_c_sum(){
	int val = 0;
	int i;
	for (i=0; i<NC; i++)
		val += c[i];

	return val;
}



/* 
==================================
 
			FILE IO		 

================================== 
*/


void packing::read_spheres(string& str){
	cout << "Reading data in from " << str << endl;

	// local vars
	int i,j,d;

	ifstream obj(str.c_str());

	if (!obj.is_open()){
		cout << "file did not open!" << endl;
		throw "file not opened!\n";
	}

	L = new double[NDIM];

	// read header info
	obj >> N;
	for (d=0; d<NDIM; d++)
		obj >> L[d];
	obj >> phi;

	cout << "N = " << N << endl;
	for (d=0; d<NDIM; d++)
		cout << "L[" << d << "] = " << L[d] << endl;
	cout << "phi = " << phi << endl;

	// read trash header
	string header(" ");
	getline(obj,header);
	cout << endl << "header #1: " << header << endl;
	getline(obj,header);
	cout << endl << "header #2: " << header << endl;
	getline(obj,header);
	cout << endl << "header #3: " << header << endl;

	// allocate memory for x,r,m
	x = new double*[N];
	r = new double[N];
	m = new double[N];

	// read data
	for (i=0; i<N; i++){
		x[i] = new double[NDIM];

		obj >> j;
		if (j != i){
			cout << "read_spheres(): data index does not match particle index, ending program." << endl;
			throw "Index mismatch!\n";
		}
		obj >> r[i];
		for (d=0; d<NDIM; d++)
			obj >> x[i][d];
		obj >> m[i];
	}
}














