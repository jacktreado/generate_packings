#ifndef PACKING_H
#define PACKING_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

class packing{
protected:
	// system information
	int N; 							// number of particles
	int NDIM; 						// dimension of space
	int DOF;						// degree of freedom per particle
	int NC;							// number of possible interparticle contacts
	double dt;						// time step
	double* L;						// array of box lengths (rectangular)	
	int seed;						// initial seed

	// particle information
	double** x;						// particle positions
	double** v;						// particle velocities
	double** F;						// forces on particles
	double** aold;					// particle acceleration from previous time step (for velocity verlet)
	double* r;						// particle radii
	double* m; 						// particle masses
	double* c;						// contact matrix (stored in 1d array)
	int* pc; 						// particle contact list

	// packing information
	double phi;						// packing fraction
	double ep; 						// energy scale for forces
	double U;						// total potential energy
	double K;						// total kinetic energy	
	int nr;							// number of rattler particles

	// neighbor list/cell list information
	double rcut; 					// multiple of radius that you check for neighbors
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
	packing(int n, int ndim, double alpha, double phi0, int seed);					// N spheres, no neighbor list
	packing(int n, int ndim, double alpha, double phi0, 
		int nc, int nnu, double rc, int seed);										// N spheres w/ neighbor list
	~packing();	

	// initialization
	void initialize_NC(){NC = (N*(N-1))/2;}
	void initialize_box(double val);
	void initialize_particles();
	void initialize_nlcl();
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
	double get_rcut() {return rcut;};
	double get_mean_vel();
	double get_com();
	double get_mean_mass();
	double get_mean_rad();
	double get_c_sum();


	// setters
	void reset_c(){ int i; for (i=0; i<NC; i++) c[i]=0; }
	void set_alpha0(double val) {alpha0 = val;};
	void set_alpha(double val) {alpha = val;};
	void set_finc(double val) {finc = val;};
	void set_fdec(double val) {fdec = val;};
	void set_falpha(double val) {falpha = val;};
	void set_dtmax(double val) {dtmax = val*dt;};
	void set_plotskip(int val) {plotskip = val;};
	void set_plotit(int val) {plotit = val;};
	void set_dt(double val) {dt = val;};
	void set_ep(double val) {ep = val;};
	void set_rcut(double val) {rcut = val;};
	void set_phi(double val) {
		int i;
		double g = pow( val/phi, (1/(double)NDIM) );
		for (i=0; i<N; i++){
			r[i] *= g;
			m[i] *= pow(g,NDIM);
		}
	};
	void set_md_time(double dt0);
	void set_rand_c(double p);
	void update_cell_g(){int i; for (i=0; i<NDIM; i++) g[i]=L[i]/NCL; };

	// file openers
	void open_xyz(std::string str){
		xyzobj.open(str.c_str());
		if (!xyzobj.is_open()){
			std::cout << "ERROR: xyzobj could not open..." << std::endl;
			std::cout << "ERROR: file str = " << str;
			throw "obj not open";
		}
	};
	void open_config(std::string str){
		configobj.open(str.c_str());
		if (!configobj.is_open()){
			std::cout << "ERROR: configobj could not open..." << std::endl;
			std::cout << "ERROR: file str = " << str;
			throw "obj not open";
		}
	};
	void open_stat(std::string str){
		statobj.open(str.c_str());
		if (!statobj.is_open()){
			std::cout << "ERROR: statobj could not open..." << std::endl;
			std::cout << "ERROR: file str = " << str;
			throw "obj not open";
		}
	};
	void open_en(std::string str){
		enobj.open(str.c_str());
		if (!enobj.is_open()){
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
	void print_cell();
	void print_neighborlist();
	void print_cell_neighbors();
	void print_nl_xyz();	
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
