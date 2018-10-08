#ifndef PACKING_H
#define PACKING_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

/// This is the packing base class. The rigidbody class inherits from it.
/**  This class contains embodies a packing experiement. */

class packing {
protected:
	// system information
	int N; 							//!< The number of particles in the experiment.
	int NDIM; 						//!< The physical dimension of the space. This can be an integer value of 1, 2, or 3
	int DOF;						//!< The number of degrees of freedom per particle
	int NC;							//!< The number of possible interparticle contacts
	double dt;						//!< The length of the time step
	double* L;						//!< Pointer to array of box lengths (rectangular)
	int seed;						//!< The initial seed

	// particle information
	double** x;						//!< Pointer to two-dimensional array containing particle positions
	double** v;						//!< Pointer to two-dimensional array containing particle velocties
	double** F;						//!< Pointer to two-dimensional array containing forces on particles
	double** aold;					//!< Pointer to two-dimensional array containing particle acceleration from previous time step (for velocity verlet)
	double* r;						//!< Particle radii
	double* m; 						//!< Particle masses
	double* c;						//!< Pointer to contact matrix (stored in 1d array)
	int* pc; 						//!< Pointer to particle contact list

	// packing information
	double phi;						//!< Packing fraction. Values range between 0 and 1
	double ep; 						//!< The energy scale for forces
	double U;						//!< The total potential energy
	double K;						//!< The total kinetic energy
	int nr;							//!< The number of rattler particles

	// neighbor list/cell list information
	double* rcut; 					// multiple of radius that you check for neighbors
	int NCL;						// number of cells along one box length
	int NCELLS;						// number of cells
	int NCN;						// number of cell neighbors (FEV neighbors on cubic lattice)
	int nnupdate;					// number of time steps to skip before updating neighbor list
	double* g;						// array of cell dimensions
	int* celln;						// number of particles in each cell
	int** cellneighbors;			// neighboring cells
	int* clabel;					// cell label for each particle
	double** cellpos;				// cell center positions
	std::vector<int>* cell;			// array of vectors: indices of atoms in each cell m
	std::vector<int>* neighborlist; // array of vectors: indices of atoms in neighbor list of atom i

	// FIRE information
	double dtmax;					// maximum time step
	double alpha;					// cg scale
	double alpha0;					// starting cg scale
	double finc;					// increase dt factor
	double fdec;					// decrease dt factor
	double falpha;					// increase alpha factor
	int np;							// number of positive P in a row

	// output information
	int plotskip;
	int plotit;
	int isjammed;
	std::ofstream xyzobj;
	std::ofstream configobj;
	std::ofstream statobj;
	std::ofstream enobj;
public:
	// constructors & destructors
	packing(int n, int dof, int nc, int s);											// N spheres, use by rigidbody class
	packing(std::string &str, int ndim, int s);										// N spheres, read from file
	packing(int n, int ndim, double alpha, double phi0,
	        int nc, int nnu, int seed);												// N spheres w/ neighbor list
	~packing();

	// initialization
	void initialize_NC() {NC = (N * (N - 1)) / 2;}
	void initialize_box(double val);
	void initialize_particles();
	void initialize_particles(int seed, double alpha, double rad);
	void initialize_nlcl();
	void nlcl_null();
	void setup_nlcl();
	void freeze_nlcl();	
	void setup_std_FIRE();
	void setup_std_sim();
	void mc_pos_init();
	void mxwl_vel_init(double t);
	void rand_vel_init(double tmp0);

	// file IO
	void read_spheres(std::string &str);

	// getters
	int get_N() {return N;};
	int get_NDIM() {return NDIM;};
	int get_DOF() {return DOF;};
	int get_seed() {return seed;};
	int get_isjammed() {return isjammed;}
	double get_dt() {return dt;};
	double get_phi() {return phi;};
	double get_ep() {return ep;};
	double get_U() {return U;};
	double get_K() {return K;};
	double get_alpha() {return alpha;};
	double get_alpha0() {return alpha0;};
	double get_finc() {return finc;};
	double get_fdec() {return fdec;};
	double get_falpha() {return falpha;};
	double get_dtmax() {return dtmax;};
	double get_distance(int p1, int p2);
	double get_distance(int p1, int p2, double xij[]);
	double get_mean_vel();
	double get_com();
	double get_mean_mass();
	double get_mean_rad();
	double get_c_sum();


	// setters
	void reset_c() { int i; for (i = 0; i < NC; i++) c[i] = 0; }
	void set_alpha0(double val) {alpha0 = val;};
	void set_alpha(double val) {alpha = val;};
	void set_finc(double val) {finc = val;};
	void set_fdec(double val) {fdec = val;};
	void set_falpha(double val) {falpha = val;};
	void set_dtmax(double val) {dtmax = val * dt;};
	void set_plotskip(int val) {plotskip = val;};
	void set_plotit(int val) {plotit = val;};
	void set_dt(double val) {dt = val;};
	void set_ep(double val) {ep = val;};
	void set_phi(double val) {
		int i;
		double g = pow( val / phi, (1 / (double)NDIM) );
		for (i = 0; i < N; i++) {
			r[i] *= g;
			m[i] *= pow(g, NDIM);
		}
	};
	void set_md_time(double dt0);
	void set_rand_c(double p);
	void update_cell_g() {int i; for (i = 0; i < NDIM; i++) g[i] = L[i] / NCL; };

	// file openers
	void open_xyz(std::string str) {
		xyzobj.open(str.c_str());
		if (!xyzobj.is_open()) {
			std::cout << "ERROR: xyzobj could not open..." << std::endl;
			std::cout << "ERROR: file str = " << str;
			throw "obj not open";
		}
	};
	void open_config(std::string str) {
		configobj.open(str.c_str());
		if (!configobj.is_open()) {
			std::cout << "ERROR: configobj could not open..." << std::endl;
			std::cout << "ERROR: file str = " << str;
			throw "obj not open";
		}
	};
	void open_stat(std::string str) {
		statobj.open(str.c_str());
		if (!statobj.is_open()) {
			std::cout << "ERROR: statobj could not open..." << std::endl;
			std::cout << "ERROR: file str = " << str;
			throw "obj not open";
		}
	};
	void open_en(std::string str) {
		enobj.open(str.c_str());
		if (!enobj.is_open()) {
			std::cout << "ERROR: enobj could not open..." << std::endl;
			std::cout << "ERROR: file str = " << str;
			throw "obj not open";
		}
	};

	// MD
	void pos_update();
	double hs(double sij, double xij);
	void hs_force();
	void hs_force_nn();
	void vel_update();
	int rmv_rattlers(int& krcrs);
	void catch_nans(int t);

	// time evolution
	void md_evol_hs(double tmp0, int NT);
	void jamming_finder(double tend, double dphi, double Utol, double Ktol);
	void jamming_finder_nn(double tend, double dphi, double Utol, double Ktol);
	void fire();
	void root_search(double& phiH, double& phiL, int& check_ratterls, int epconst,
	                 int nr, double dphi0, double Ktol, double Utol, int t);

	// neighbor list/cell list
	void init_rcut();
	void scale_rcut(double s);
	void update_nlcl(int t);
	void print_cell();
	void print_neighborlist();
	void print_cell_neighbors();
	void print_nl_xyz();
	void print_nl_xyz(int p1, int p2);
	void print_all_nl_xyz();
	void print_cell_pos();
	void print_clabel();
	void print_celln();
	void add_cell(int m, int i);
	void reset_cells();
	void add_neighbor(int i, int j);
	void reset_neighborlist(int i);
	void update_celln();
	void get_cell_neighbors();
	void get_cell_positions();
	int get_new_cell(int i);
	void update_cell();
	void update_neighborlist();
	void nlcl_check(int t);

	// growth
	void scale_sys(double dphi);

	// measurements
	void calc_vacf(int NT, int vsave, std::ofstream& obj);
	void get_vacf(int t, int vsave, std::vector<double>** vlist, std::vector<double>& numer, std::vector<double>& denom);
	void finish_vacf(int NT, int vsave, std::vector<double>& numer, std::vector<double>& denom, std::vector<double>& vacf);
	void print_vacf(int NT, int vsave, std::vector<double>& vacf, std::ofstream& obj);

	void calc_dm(std::string& fstr);
	double perturb_single_particle(int i, int k, double h);
	double perturb_two_particles(int i, int j, int k, int l, double h);

	// printers
	void md_monitor(int t, int nr, double phiH, double phiL);
	void monitor_scale(double dphi, double phiH, double phiL);
	void print_xyz();
	void print_c_mat();
	void print_c_mat(std::ofstream &obj);
	void print_c_data(std::ofstream &obj);
	void print_pc(std::ofstream &obj, int w);
	void print_vars();
	void print_config();
	void print_stat();
	void monitor_header(int t);

};
#endif
