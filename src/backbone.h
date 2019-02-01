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
	// backbone atom ids
	int nid;
	int caid;
	int cid;

	// spring constants
	double kbl;	
	double kba;
	double kda;

	// angle values
	double* theta;
	double* eta;
	double* phi_da;
	double* psi_da;
	double* omega_da; 

	// rest lengths/angles
	double* l0; 
	double* theta0;
	double* eta0;
	double* phi0_da;
	double* psi0_da;
	double* omega0_da;
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

	double get_theta(int r);
	double get_eta(int r);
	double get_phi_da(int r);
	double get_psi_da(int r);
	double get_omega_da(int r);

	// simulation
	void bb_free_md(double T0, int NT, int nnu);

	// topology relax
	void top_relax();	

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
};
#endif