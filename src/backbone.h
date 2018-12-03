#ifndef BACKBONE_H
#define BACKBONE_H

#include "Quaternion.h"
#include "packing.h"
#include "rigidbody.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

class backbone : public rigidbody {
protected:
	// connection list
	int** cnx;

	// spring constants
	double kbl;	
	double kba;
	double kda;
public:
	// constructor and destructor	
	backbone(string &bbstr, int n, int dof, int nc, int s);
	~backbone();

	// initialization
	void initialize_cnx(string &bbstr, int NLINES);

	// setters/getters
	void set_kbl (double k) {kbl = k;};
	void set_kba (double k) {kba = k;};
	void set_kda (double k) {kda = k;};

	// topology relax
	void top_relax();

	// force update
	double bb_force_update();
};

// CONSTRUCTOR

backbone::backbone(string &bbstr, int n, int dof, int nc, int s) : rigidbody(string &bbstr, int n, int dof, int nc, int s) {
	cout << "in backbone constructor..." << endl;

	// initialize connection list cnx
	int NATOT = this->get_Natot();
	this->initialize_cnx(bbstr,NATOT);

	cout << "connections established, backbone constructor complete." << endl;
}

// DESTRUCTOR

backbone::~backbone(){
	// local variables
	int i;

	// free memory for cnx
	for (i=0; i<N; i++)
		delete [] cnx[i];
	delete [] cnx;
}




// INITIALIZATION

void backbone::initialize_cnx(string &bbstr, int NATOT){
	// local variables
	int i;
	ifstream obj(bbstr.c_str());

	// allocate memory
	cnx = new int*[N];
	cnx[0] = new int[1];
	for (i=1; i<N-1; i++)
		cnx[i] = new int[2];
	cnx[N-1] = new int[1];


	// loop over lines in input file, skip to backbone info
	string trash;
	NLINES = 6+NATOT;
	for (i=0; i<NLINES; i++) { getl(obj,trash); }

	// loop over residues, populate connections
	obj >> cnx[0][0];	
	for (i=1; i<N-1; i++)
		obj >> cnx[i][0] >> cnx[i][1];
	obj >> cnx[N-1][0];
}


// BACKBONE: RELAX TOPOLOGY

void backbone::top_relax(){
	// local variables
	int i;
	double Ubb = 0;
	double Utol = 1e-16;

	// loop while backbone has a high potential energy
	while (Ubb > Utol){
		// advance quaternions, positions
		this->verlet_first();

		// update backbone forces
		Ubb = this->bb_force_update();

		// include FIRE relaxation
		this->rb_fire();

		// advance angular momentum
		this->verlet_second();
	}

}


// BACKBONE FORCE UPDATE

double backbone::bb_force_update(){
	// local variables
	int i,ip1;
	double U = 0; 

	// loop over particles, calculate forces due to backbone connections
	for (i=0; i<N-1; i++){
		// get index of next residue
		ip1 = i+1;

		// bond length force
		U = U + this->bl_force(i,ip1);

		// bond angle force
		U = U + this->ba_force(i,ip1);

		// dihedral angle force
		U = U + this->da_force(i,ip1);
	}

	return U;
}



// bond length force
double backbone::bl_force(int i, int j){
	// local variables
	double dx,dy,dz;

	return U;
}


// bond angle force
double backbone::bl_force(int i, int j){
	// local variables
	double dx,dy,dz;

	return U;
}


// dihedral angle' force
double backbone::bl_force(int i, int j){
	// local variables
	double dx,dy,dz;

	return U;
}





#endif