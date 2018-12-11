#include "Quaternion.h"
#include "packing.h"
#include "rigidbody.h"
#include "backbone.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

const double PI = atan(1)*4;

// CONSTRUCTOR

backbone::backbone(string &bbstr, int n, int dof, int nc, int s) : rigidbody(bbstr, n, dof, nc, s) {
	cout << "in backbone constructor..." << endl;

	// initialize connection list cnx
	int NATOT = this->get_Natot();
	this->initialize_bb(bbstr,NATOT);

	// set spring constants to 1
	kbl = 1.0;
	kba = 1.0;
	kda = 1.0;

	// set rest lengths to 100% of contact distance
	l0 = 2;
	theta0 = 0.2*PI;
	da0 = 0.2*PI;

	cout << "connections established, backbone constructor complete." << endl;
}

// DESTRUCTOR

backbone::~backbone(){
	// local variables
	int i;

	// free memory for cnx
	for (i=0; i<N; i++){
		delete [] cnx[i];
		delete [] bbid[i];
	}
	delete [] cnx;
	delete [] bbid;
}




// INITIALIZATION

void backbone::initialize_bb(string &bbstr, int NATOT){
	// local variables
	int i,NLINES;
	ifstream obj(bbstr.c_str());

	// allocate memory
	cnx = new int*[N];
	bbid = new int*[N];
	cnx[0] = new int[1];
	bbid[0] = new int[2];
	for (i=1; i<N-1; i++){
		cnx[i] = new int[2];
		bbid[i] = new int[2];
	}
	cnx[N-1] = new int[1];
	bbid[N-1] = new int[2];

	// loop over lines in input file, skip to backbone info
	string trash;
	NLINES = 5+NATOT;
	cout << "Getting backbone info from input file now..." << endl;
	cout << "TRASH:" << endl;
	for (i=0; i<NLINES; i++) { getline(obj,trash); cout << trash << endl; }
	cout << "Now getting the good stuff..." << endl;

	// loop over residues, populate connections
	obj >> bbid[0][0] >> bbid[0][1] >> cnx[0][0];
	cout << bbid[0][0] << bbid[0][1] << cnx[0][0] << endl;
	for (i=1; i<N-1; i++){
		obj >> bbid[i][0] >> bbid[i][1] >> cnx[i][0] >> cnx[i][1];
		cout << bbid[i][0] << bbid[i][1] << cnx[i][0] << cnx[i][1] << endl;
	}
	obj >> bbid[N-1][0] >> bbid[N-1][1] >> cnx[N-1][0];
	cout << bbid[N-1][0] << bbid[N-1][1] << cnx[N-1][0] << endl;

	// close ifstream object
	obj.close();
}


// BACKBONE: RELAX TOPOLOGY

void backbone::top_relax(){
	// local variables
	int k = 0;
	int kmax = 1e3;
	double Ubb = 0;
	double Utol = 1e-16;
	double T0 = 1e-2;

	// initialize velocities
	this->rand_vel_init(T0);

	// get backbone forces
	Ubb = this->bb_force_update();

	// loop while backbone has a high potential energy
	while (Ubb > Utol && k < kmax){
		// advance quaternions, positions
		this->verlet_first();

		// update backbone forces
		Ubb = this->bb_force_update();
		U = Ubb;

		// include FIRE relaxation
		this->rb_fire();

		// advance angular momentum
		this->verlet_second();

		// print some stuff
		if (k % plotskip == 0){
			this->monitor_header(k);
			this->rigidbody_md_monitor();
			cout << endl;
			cout << endl;
		}

		// update iterate
		k++;
	}

	// report on success or failure of loop
	if (k == kmax){
		cout << "ERROR: Loop did not complete, ending..." << endl;
		throw;
	}
	else
		cout << "Loop completed! topology of backbone sufficiently relaxed..." << endl;

}


// BACKBONE FORCE UPDATE

double backbone::bb_force_update(){
	// local variables
	int i,ip1;
	double Ubb = 0; 

	// loop over particles, calculate forces due to backbone connections
	for (i=0; i<N; i++){
		// bond length force
		if (i > 0)
			Ubb = Ubb + this->bl_force(i);

		// // bond angle force
		// if (i > 1)
		// 	U = U + this->ba_force(i);

		// // dihedral angle force
		// if (i > 2)
		// 	U = U + this->da_force(i);
	}

	return Ubb;
}



// bond length force
double backbone::bl_force(int i){
	// local variables
	double dij[NDIM];
	double vfij[NDIM];
	double rij,fij,r0,Ubb;
	double qix,qiy,qiz,qjx,qjy,qjz;
	double tix,tiy,tiz,tjx,tjy,tjz;
	int im1,nid,cid,d;

	// get attached monomer
	im1 = i-1;
	cid = bbid[im1][1];
	nid = bbid[i][0];

	// calculate distance
	rij = this->get_atomic_distance(im1,i,cid,nid,dij);

	// calculate scalar force
	r0 = ar[im1][cid]+ar[i][nid];
	r0 = l0*r0;
	fij = kbl*(1-(r0/rij));

	// update vectorial force
	for (d=0; d<NDIM; d++){
		vfij[d] = fij*dij[d];
		F[im1][d] += vfij[d];
		F[i][d] -= vfij[d];
	}

	// branch vectors to atoms
	qix = xW[i][cid][0];
	qiy = xW[i][cid][1];
	qiz = xW[i][cid][2];

	qjx = xW[im1][nid][0];
	qjy = xW[im1][nid][1];
	qjz = xW[im1][nid][2];

	// update torques
	TqW[im1][0] += qiy * vfij[2] - qiz * vfij[1];
	TqW[im1][1] += -qix * vfij[2] + qiz * vfij[0];
	TqW[im1][2] += qix * vfij[1] - qiy * vfij[0];

	TqW[i][0] += -qjy * vfij[2] + qjz * vfij[1];
	TqW[i][1] += qjx * vfij[2] - qjz * vfij[0];
	TqW[i][2] += -qjx * vfij[1] + qjy * vfij[0];

	// potential
	Ubb = 0.5*kbl*pow(rij-r0,2);

	// diagnostic output
	// cout << "** Distances and forces between (" << im1 << "," << cid << ") and (" << i << "," << nid << ") = " << endl;
	// cout << "distance: ";
	// cout << setw(10) << rij/r0;
	// cout << setw(10) << dij[0];
	// cout << setw(10) << dij[1];
	// cout << setw(10) << dij[2];
	// cout << endl;
	// cout << "forces: ";
	// cout << setw(10) << vfij[0];
	// cout << setw(10) << vfij[1];
	// cout << setw(10) << vfij[2];
	// cout << endl << endl;


	return Ubb;
}


// bond angle force
double backbone::ba_force(int i){
	// local variables
	double dx,dy,dz,Ubb;

	return Ubb;
}


// dihedral angle' force
double backbone::da_force(int i){
	// local variables
	double dx,dy,dz,Ubb;

	return Ubb;
}
