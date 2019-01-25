#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include "Quaternion.h"
#include "packing.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

/**
 * The rigidbody class inherits from the packing class and contains
 * functionality to simulate a packing experiement consisting of rigidbodies.
 */

class rigidbody : public packing{
protected:
	double Krot; //!< rotational kinetic energy
	double LCON; //!< local angular momentum

	// rigid body info (quaternion, angles)
	double*** xW;							//!< atomic World frame position
	double*** xM;							//!< atomic Molecule frame position									
	Quaternion* q;							//!< main quaternion
	Quaternion* qhalf;						//!< half step quaternion
	Quaternion* qnew;						//!< update quaternion
	Quaternion* qdot;						//!< quaternion derivative
	double* eulang1;						//!< 1st euler angle: azimuth (phi)
	double* eulang2;						//!< 2nd euler angle: polar (theta)
	double* eulang3;						//!< 3rd euler angle: spin (psi)
	double** Inn;							//!< Moment of inertia tensor (diagnolized)

	// velocity-verlet rigid body dynamics		
	double** wW;								//!< world frame angular velocity vector
	double** wM;								//!< molecule frame angular velocity vector	
	double** TqW;								//!< World frame torque
	double** TqM;								//!< Molecule frame torque
	double** LW;								//!< World frame angular momentum
	double** LWhalf;							//!< half step World frame angular momentum
	double** LM;								//!< Molecule frame angular momentum
	double** LMhalf;							//!< half step Molecule frame angular momentum
	double** LMdot;								//!< Molecule frame angular momenum time derivative

	// per-atom info
	int* Na;									//!< number of atoms per particle
	int* ac;									//!< atomic contact array, counting bump-bump contacts
	int* cm;									//!< contact multiplicity matrix
	double** ar;								//!< radii of each atom	
public:
	rigidbody(string &rbstr, int n, int dof, int nc, int s);
	~rigidbody();

	// getters
	void reset_cm();
	void update_phi();
	void update_euler();
	int get_ac_sum();
	double get_atomic_distance(int i, int j, int ai, int aj, double aij[]);
	double get_Natot();
	double get_LWX();
	double get_LWY();
	double get_LWZ();

	// initialization
	void get_file_header(string &rbstr);
	void initialize_particles();
	void initialize_dynamics();
	void read_in_info(string &rbstr);
	void initialize_quaternions();

	// quaternions
	void q_step();
	void rotation_W2M();
	void euler_q();
	void pos_frot();
	void pos_brot();

	// Molecular Dynamics
	void free_md(double tmp0, int NT, int nnu);
	void free_fire(double tmp0, double Utol, double tend, int nnu);
	void single_md(int t);
	void single_fire(int t);

	// MD components
	void verlet_first();
	void verlet_second();	
	void force_update();
	int rmv_rattlers(int& krcrs);
	int id_rattlers();	

	// packing growth
	void rb_jamming_finder(double tmp0, int NT, double dphi, double Utol, double Ktol);
	void rb_jamming_precise(double tphiold, double tmp0, int NT, double Utol, double Ktol);
	void rb_jamming_easy(double tmp0, int NT, double dphi, double Utol, double Ktol);
	void get_U(double Ktol, int& nr);
	void rb_scale(double phinew);
	void rb_root_search(double& phiH, double& phiL, int& check_rattlers, int epconst, int nr, double dphi0, double Ktol, double &Utol, int t);
	void rb_easy(double& phiH, double& phiL, int& check_rattlers, int &epconst, int nr, double dphi0, double Ktol, double &Utol, int t, bool &min);
	void rb_fire();

	// Printing
	void print_stat();
	void print_ac(std::ofstream& obj, int w);
	void print_config();
	void monitor_scale(double dphi, double phiL, double phiH);
	void rigidbody_md_monitor();
	void rigidbody_xyz();
	void rigidbody_xyz(int p1, int p2);
	void rigidbody_print_vars();

};
#endif














