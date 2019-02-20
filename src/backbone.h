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
	// energy termns
	double Ubb;
	double Ubl;
	double Uba;
	double Uda;
	double Usteric;

	// backbone atom ids
	int nid;
	int caid;
	int cid;

	// spring constants
	double kbl;	
	double kba;
	double kda;

	// cosine of angle values
	double* ctheta;
	double* ceta;
	double* cphi;
	double* cpsi;
	double* comega; 

	// rest lengths/cosines
	double* l0; 
	double* ctheta0;
	double* ceta0;
	double* cphi0;
	double* cpsi0;
	double* comega0;
public:
	// constructor and destructor	
	backbone(string &bbstr, int n, int dof, int nc, int s);
	~backbone();

	// initialization
	void bl0_init();
	void bl0_init(double val);

	void ba0_init();
	void ba0_init(double val);

	void da0_init();
	void da0_init(double val);

	// setters/getters
	void set_kbl(double val) {kbl = val;};
	void set_kba(double val) {kbl = val;};
	void set_kda(double val) {kbl = val;};
	void set_angles();

	double get_ctheta(int r);
	double get_ceta(int r);
	double get_cphi(int r);
	double get_cpsi(int r);
	double get_comega(int r);

	// simulation
	void bb_free_md(double T0, int NT, int nnu);

	// topology relax
	void top_relax(int kmax);	

	// scaling
	void rb_scale(double phinew);

	// force update
	double steric_force_update();
	double bb_force_update (int choice);
	double bl_force (int i);
	double ba_force (int i);
	double da_force (int i);
	double dotp(double a[], double b[]);

	// printing
	void print_angles();
	void bb_md_monitor();
};
#endif