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

using namespace std;

class backbone : public rigidbody {
protected:
	// connection list
	int** cnx;

	// bb atom indices
	int** bbid;

	// spring constants
	double kbl;	
	double kba;
	double kda;	

	// rest lengths/angles
	double l0; 
	double theta0;
	double da0;
public:
	// constructor and destructor	
	backbone(string &bbstr, int n, int dof, int nc, int s);
	~backbone();

	// initialization
	void initialize_bb(string &bbstr, int NLINES);

	// setters/getters
	void set_l0 (double l) {l0 = l;};
	void set_kbl (double k) {kbl = k;};
	void set_kba (double k) {kba = k;};
	void set_kda (double k) {kda = k;};

	double get_l0 () {return l0;}

	// topology relax
	void top_relax ();

	// force update
	double bb_force_update ();
	double bl_force (int i);
	double ba_force (int i);
	double da_force (int i);
};
#endif