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


/****************************
*							*
*		CONSTRUCTORS		*
*	   		AND				*
*		DESTRUCTORS			*
*							*
*****************************/


// CONSTRUCTOR

// test constructor
backbone::backbone(string &bbstr, int n, int dof, int nc, int s) : rigidbody(bbstr, n, dof, nc, s) {
	cout << "in test backbone constructor..." << endl;
	cout << endl;

	// local variables
	int i;
	double len0=0, th0=0, da0=0;

	// set indices
	nid = 0;
	caid = 1;
	cid = 2;

	// set spring constants to 1
	kbl = 1.0;
	kba = 1.0;
	kda = 1.0;

	// allocate memory
	cout << "allocating memory to angle and rest angle arrays..." << endl;

	// arrays of values
	theta = new double [N];
	eta = new double [N];
	phi_da = new double [N];
	psi_da = new double [N];
	omega_da = new double [N];

	// set initial values to -1
	for (i=0; i<N; i++){
		theta[i] = -1.0;
		eta[i] = -1.0;
		phi_da[i] = -1.0;
		psi_da[i] = -1.0;
		omega_da[i] = -1.0;
	}

	// rest angles
	l0 = new double [N];
	theta0 = new double [N];
	eta0 = new double [N];
	phi0_da = new double [N];
	psi0_da = new double [N];
	omega0_da = new double [N];

	// set rest lengths/angle
	cout << "initializing values of angle arrays to std values; ";	
	len0 = 1.0;					// fraction of contact distance
	this->bl0_init(len0);
	th0 = 0.6*PI;				// theta0 (everything)
	this->ba0_init(th0);	
	da0 = 0.6*PI;				// da0 (everything)
	this->da0_init(da0);
	cout << "bl0 = " << len0 << "; ";
	cout << "ba0 = " << th0 << "; ";
	cout << "da0 = " << da0 << endl;
	cout << "PI = " << PI << endl;

	cout << endl;
	cout << "test backbone constructor complete." << endl;
	cout << endl << endl;
}

// DESTRUCTOR
backbone::~backbone(){
	// free memory for values
	delete [] theta;
	delete [] eta;
	delete [] phi_da;
	delete [] psi_da;
	delete [] omega_da;

	// free memory for rest values
	delete [] theta0;
	delete [] eta0;
	delete [] phi0_da;
	delete [] psi0_da;
	delete [] omega0_da;
}


/****************************
*							*
*							*
*	   INITIALIZATION		*
*							*
*							*
*****************************/



void backbone::bl0_init(double val){
	int i;

	// loop over indices
	for (i=0; i<N; i++)
		l0[i] = val;
}

void backbone::ba0_init(double val){
	int i;

	// loop over indices
	for (i=0; i<N; i++){
		theta0[i] = val;
		eta0[i] = val;
	}
}

void backbone::da0_init(double val){
	int i;

	// loop over indices
	for (i=0; i<N; i++){
		phi0_da[i] = val;
		psi0_da[i] = val;
		omega_da[i] = val;
	}
}












/****************************
*							*
*							*
*	SETTERS AND GETTERS		*
*							*
*							*
*****************************/


// SETTERS

// update angles based on particle positions
void backbone::set_angles(){
	// local variables
	int i;

	// loop over particles, measure angles
	eta[0] = this->get_eta(0);					// eta depends on i+1, so has i=0 value
	psi_da[0] = this->get_psi_da(0);			// psi depends on i+1, so has i=0 value
	for (i=1; i<N-1; i++){
		eta[i] = this->get_eta(i);
		theta[i] = this->get_theta(i);		
		phi_da[i] = this->get_phi_da(i);
		psi_da[i] = this->get_psi_da(i);
		omega_da[i] = this->get_omega_da(i);
	}
	theta[N-1] = this->get_theta(N-1);			// theta depends on i-1, so has i=N-1 value
	phi_da[N-1] = this->get_phi_da(N-1);		// phi depends on i-1, so has i=N-1 value
	omega_da[N-1] = this->get_omega_da(N-1);	// omega depends on i-1, so has i=N-1 value
}



// GETTERS

// get theta angle for residue r
double backbone::get_theta(int r){
	// connection vectors
	double v1[NDIM];
	double v2[NDIM];
	double num,v1norm,v2norm,denom,costheta,val;
	int rm1,d;

	// set indices
	rm1 = r-1;

	// get vectors and norms
	v1norm = this->get_atomic_distance(rm1,r,cid,nid,v1);
	v2norm = this->get_atomic_distance(r,r,nid,caid,v2);

	// numerator (use definition of C)
	num = this->dotp(v1,v2);

	// denominator
	denom = v1norm*v2norm;

	// get cosine
	costheta = -1*num/denom;

	// use inverse cosine to get eta
	val = acos(costheta);

	// return value
	return val;
}

// get eta angle for residue r
double backbone::get_eta(int r){
	// connection vector
	double v1[NDIM];
	double v2[NDIM];
	double num,v1norm,v2norm,denom,coseta,val;
	int rp1,d;

	// set indices
	rp1 = r+1;

	// get vectors and norms
	v1norm = this->get_atomic_distance(r,r,caid,cid,v1);
	v2norm = this->get_atomic_distance(r,rp1,cid,nid,v2);

	// numerator (use definition of C)
	num = this->dotp(v1,v2);

	// denominator
	denom = v1norm*v2norm;

	// get cosine
	coseta = -1*num/denom;

	// use inverse cosine to get eta
	val = acos(coseta);

	// return value
	return val;
}

// get phi dihedral angle for residue r (WORK IN PROGRESS)
double backbone::get_phi_da(int r){
	return 0;
}

// get psi dihedral angle for residue r (WORK IN PROGRESS)
double backbone::get_psi_da(int r){
	return 0;
}

// get omega dihedral angle for residue r (WORK IN PROGRESS)
double backbone::get_omega_da(int r){
	return 0;
}



/****************************
*							*
*							*
*		 SIMULATION		    *
*							*
*							*
*****************************/



// BACKBONE: RELAX TOPOLOGY
void backbone::top_relax(){
	// local variables
	int k = 0;
	int kmax = 2e3;
	double Ubb = 0;
	double Usteric = 0;
	double Utol = 1e-16;
	double T0 = 1e-2;

	// initialize velocities
	this->rand_vel_init(T0);

	// get backbone forces
	this->set_angles();
	this->print_angles();
	Ubb = this->bb_force_update(0);
	U = Ubb;

	// loop while backbone has a high potential energy
	while (k < kmax){
		// advance quaternions, positions
		this->verlet_first();

		// update angles
		this->set_angles();

		// update backbone forces
		Ubb = this->bb_force_update(1);
		U = Ubb;

		// include FIRE relaxation
		this->rb_fire();

		// advance angular momentum
		this->verlet_second();

		// print some stuff
		if (k % plotskip == 0){
			this->monitor_header(k);
			this->rigidbody_md_monitor();
			this->print_angles();
			cout << endl;
			cout << endl;
		}

		// update iterate
		k++;
	}

	cout << "Loop completed! topology of backbone sufficiently relaxed..." << endl;

	// Now do hard sphere steric relaxation with bond length constraints
	k = 0;	
	do {
		// advance quaternions, positions
		this->verlet_first();

		// update angles
		this->set_angles();

		// update backbone forces
		Usteric = 0;
		Ubb = 0;
		Usteric = this->steric_force_update();
		// Ubb = this->bb_force_update(0);
		U = Ubb + Usteric;

		// include FIRE relaxation
		this->rb_fire();

		// advance angular momentum
		this->verlet_second();

		// print some stuff
		if (k % plotskip == 0){
			this->monitor_header(k);
			this->rigidbody_md_monitor();
			this->print_angles();
			cout << endl;
			cout << endl;
		}

		// update iterate
		k++;
	} while (U > Utol && k < kmax);

	// report on success or failure of loop
	if (k == kmax){
		cout << "ERROR: Loop did not complete, ending..." << endl;
		throw;
	}
	else
		cout << "Loop completed! topology of backbone sufficiently relaxed..." << endl;
}

void backbone::bb_free_md(double tmp0, int NT, int nnu) {
	int t;
	double Ubb,Usteric;

	// initialize velocities
	this->rand_vel_init(tmp0);

	// update nearest neighbor update
	nnupdate = nnu;

	cout << "===== BEGINNING TIME LOOP ======" << endl;
	cout << "NT = " << NT << endl;
	cout << "dt = " << dt << endl;

	for (t = 0; t < NT; t++) {
		// if NLCL, update
		if (NCL > -1 && t % nnupdate == 0) {
			this->update_cell();
			this->update_neighborlist();
		}

		// advance quaternions, positions
		this->verlet_first();

		// update angles
		this->set_angles();

		// update forces
		Usteric = 0;
		Ubb = 0;

		Usteric = this->steric_force_update();
		Ubb = this->bb_force_update(1);
		U = Usteric + Ubb;		


		// advance angular momentum
		this->verlet_second();

		// output information
		if (t % plotskip == 0) {
			this->monitor_header(t);
			this->rigidbody_md_monitor();
			this->print_angles();
			cout << endl;
			cout << endl;
		}
	}
}

void backbone::rb_scale(double phinew){ 
	int i; 
	double s,invdim;

	// get scaling parameter
	invdim = (double)(1/NDIM);
	s = pow(phinew/phi, invdim);

	// scale rest lengths
	for(i=0; i<N-1; i++) {
		l0[i] *= s;
	}

	// call rigibody scale
	this->rigidbody::rb_scale(phinew); 
}






/****************************
*							*
*							*
*		  FORCES &		    *
*		  TORQUES			*
*							*
*							*
*****************************/


// steric force, no overlap force between Ci-1 and Ni (for i = 1,...,N-1)
double backbone::steric_force_update(){
	int i, j, jj, ai, aj, d, cind, M, pc_found;
	double sij, rij, dx, da;
	double qix, qiy, qiz, qjx, qjy, qjz;
	double tix, tiy, tiz, tjx, tjy, tjz;
	double fix, fiy, fiz;
	double aij[NDIM];
	double fij[NDIM];

	// reset Usteric, LCON
	double Usteric = 0;
	this->reset_c();
	this->reset_cm();

	// implement NLCL if there are positive number of cells
	if (NCL > 0) {
		for (i = 0; i < N; i++) {
			M = neighborlist[i].size();
			for (jj = 0; jj < M; jj++) {
				tix = 0; tiy = 0; tiz = 0;
				tjx = 0; tjy = 0; tjz = 0;
				fix = 0; fiy = 0; fiz = 0;

				// new particle contact not yet found
				pc_found = 0;

				// get neighbor from neighborlist
				j = neighborlist[i].at(jj);	

				// contact matrix index
				cind = N * i + j - ((i + 1) * (i + 2)) / 2;

				// get contact distance sij
				sij = r[j] + r[i];

				// get distance between particles
				dx = this->get_distance(i, j);

				if (dx < sij) {
					for (ai = 0; ai < Na[i]; ai++) {
						for (aj = 0; aj < Na[j]; aj++) {

							// check if overlapping Ci-1 and Ni
							if (j-i == 1 && ai == cid && aj == nid)
								continue;
							else if (i-j == 1 && ai == nid && aj == cid)
								continue;

							// contact distance
							rij = ar[i][ai] + ar[j][aj];

							// get distance between atoms
							da = this->get_atomic_distance(i, j, ai, aj, aij);

							// if true, atoms are overlapping, so calc force and torquez
							if (da < rij) {
								// update contact forces
								for (d = 0; d < NDIM; d++) {
									fij[d] = (this->hs(rij, da)) * aij[d];
									F[i][d] += fij[d];
									F[j][d] -= fij[d];
								}

								// branch vectors to atoms
								qix = xW[i][ai][0];
								qiy = xW[i][ai][1];
								qiz = xW[i][ai][2];

								qjx = xW[j][aj][0];
								qjy = xW[j][aj][1];
								qjz = xW[j][aj][2];

								// update torques
								tix += qiy * fij[2] - qiz * fij[1];
								tiy += -qix * fij[2] + qiz * fij[0];
								tiz += qix * fij[1] - qiy * fij[0];

								tjx += -qjy * fij[2] + qjz * fij[1];
								tjy += qjx * fij[2] - qjz * fij[0];
								tjz += -qjx * fij[1] + qjy * fij[0];

								// update net force due to j
								fix += fij[0];
								fiy += fij[1];
								fiz += fij[2];

								// update contact list
								if (pc_found == 0) {
									pc[i]++;
									pc[j]++;
									c[cind] = 1;
									cm[cind] = 1;
									pc_found = 1;
								}
								else
									cm[cind]++;
								ac[i]++;
								ac[j]++;								

								// update potential energy
								Usteric += (ep / 2) * pow(1 - da / rij, 2);
							}
						}
					}

					// update torque on i & j
					TqW[i][0] += tix;
					TqW[i][1] += tiy;
					TqW[i][2] += tiz;

					TqW[j][0] += tjx;
					TqW[j][1] += tjy;
					TqW[j][2] += tjz;
				}
			}
		}
	}

	// else, just do N(N-1)/2 loop when checking for contacts
	else {
		for (i = 0; i < N; i++) {
			for (j = i + 1; j < N; j++) {
				tix = 0; tiy = 0; tiz = 0;
				tjx = 0; tjy = 0; tjz = 0;
				fix = 0; fiy = 0; fiz = 0;

				// new particle contact not yet found
				pc_found = 0;

				// contact matrix index
				cind = N * i + j - ((i + 1) * (i + 2)) / 2;

				// get contact distance sij
				sij = r[j] + r[i];

				// get distance between particles
				dx = this->get_distance(i, j);

				// if true, residues close by, check atomic overlaps
				if (dx < sij) {
					for (ai = 0; ai < Na[i]; ai++) {
						for (aj = 0; aj < Na[j]; aj++) {

							// check if overlapping Ci-1 and Ni
							if (j-i == 1 && ai == cid && aj == nid)
								continue;

							// contact distance
							rij = ar[i][ai] + ar[j][aj];

							// get distance between atoms
							da = this->get_atomic_distance(i, j, ai, aj, aij);

							// if true, atoms are overlapping, so calc force and torquez
							if (da < rij) {
								// update contact forces
								for (d = 0; d < NDIM; d++) {
									fij[d] = (this->hs(rij, da)) * aij[d];
									F[i][d] += fij[d];
									F[j][d] -= fij[d];
								}								

								// update contact list
								if (pc_found == 0) {
									pc[i]++;
									pc[j]++;
									c[cind] = 1;
									cm[cind] = 1;
									pc_found = 1;
								}
								else
									cm[cind]++;
								ac[i]++;
								ac[j]++;

								// branch vectors to atoms
								qix = xW[i][ai][0];
								qiy = xW[i][ai][1];
								qiz = xW[i][ai][2];

								qjx = xW[j][aj][0];
								qjy = xW[j][aj][1];
								qjz = xW[j][aj][2];

								// update torques
								tix += qiy * fij[2] - qiz * fij[1];
								tiy += -qix * fij[2] + qiz * fij[0];
								tiz += qix * fij[1] - qiy * fij[0];

								tjx += -qjy * fij[2] + qjz * fij[1];
								tjy += qjx * fij[2] - qjz * fij[0];
								tjz += -qjx * fij[1] + qjy * fij[0];

								// update net force due to j
								fix += fij[0];
								fiy += fij[1];
								fiz += fij[2];		

								// update potential energy
								Usteric += (ep / 2) * pow(1 - da / rij, 2);
							}							
						}
					}

					// update torque on i & j
					TqW[i][0] += tix;
					TqW[i][1] += tiy;
					TqW[i][2] += tiz;

					TqW[j][0] += tjx;
					TqW[j][1] += tjy;
					TqW[j][2] += tjz;
				}
			}
		}
	}
	return Usteric;
}


// BACKBONE FORCE UPDATE

// choice can choose which interaction to turn on:
// choice = 0 -> only BL force
// choice = 1 -> only BA force
// choice = 2 -> only DA force
double backbone::bb_force_update(int choice){
	// local variables
	int i,ip1;
	double Ubb = 0;

	// loop over particles, calculate forces due to backbone connections
	for (i=0; i<N; i++){
		// bond length force (done pairwise)
		if (i > 0)
			Ubb = Ubb + this->bl_force(i);

		// bond angle force
		if (choice > 0)		
			Ubb = Ubb + this->ba_force(i);

		// // dihedral angle force
		if (choice > 1)
			Ubb = Ubb + this->da_force(i);
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
	int im1,d;

	// get attached monomer
	im1 = i-1;

	// calculate distance
	rij = this->get_atomic_distance(im1,i,cid,nid,dij);

	// calculate scalar force
	r0 = ar[im1][cid]+ar[i][nid];
	r0 = l0[im1]*r0;
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

	return Ubb;
}


// bond angle force

// forces determined by angle i, i+1, AND i+2
double backbone::ba_force(int i){
	// local variables
	int im1,ip1;						// particle and atomic indices
	int d;								// index for dimensions
	double Ubb,Utheta,Ueta;				// potential energies
	double fN[NDIM];					// force vector on N atom
	double fC[NDIM];					// force vector on C atom
	double qNx,qNy,qNz,qCx,qCy,qCz;		// branch vectors to atomic coordinates

	// connection vectors
	double wim1[NDIM];
	double ui[NDIM];
	double vi[NDIM];
	double wi[NDIM];
	double uip1[NDIM];
	double vip1[NDIM];

	// magnitude of connection vectors
	double wim1_mag, ui_mag, vi_mag, wi_mag, uip1_mag, vip1_mag;

	// dot product terms
	double Cww_im1, Cwu_im1, Cuu_i, Cvv_i, Cuv_i, Cww_i, Cwu_i, Cuu_ip1, Cvv_ip1, Cuv_ip1;

	// composite terms
	double Dwu_im1, Duv_i, Dwu_i, Duv_ip1;

	// constants in force calculation
	double K1, K2, K3, K4;	// force Konstants
	double R1, R2, R3, R4;	// Ratio of constants

	// get residue indices
	im1 = i-1;
	ip1 = i+1;	

	// calculate connection vectors
	
	// vectors with prevous atomic positions (do not exist for i = 0)
	if (i>0){
		// wim1[d] = xW[im1][cid][d] - xW[im1][caid][d];
		wim1_mag = this->get_atomic_distance(im1,im1,caid,cid,wim1);
		// ui[d] = xW[i][nid][d] - xW[im1][cid][d];
		ui_mag = this->get_atomic_distance(im1,i,cid,nid,ui);
	}

	// vectors with atoms of residue i
	// vi[d] = xW[i][caid][d] - xW[i][nid][d];
	vi_mag = this->get_atomic_distance(i,i,nid,caid,vi);
	// wi[d] = xW[i][cid][d] - xW[i][caid][d];
	wi_mag = this->get_atomic_distance(i,i,caid,cid,wi);

	// vectors with next atom positions (do not exist for i = N-1)
	if (i<N-1){
		// uip1[d] = xW[ip1][nid][d] - xW[i][cid][d];
		uip1_mag = this->get_atomic_distance(i,ip1,cid,nid,uip1);
		// vip1[d] = xW[ip1][caid][d] - xW[ip1][nid][d];
		vip1_mag = this->get_atomic_distance(ip1,ip1,nid,caid,vip1);
	}

	// calculate dot product terms
	// involving i-1
	if (i>0){
		// Cww_im1 = this->dotp(wim1,wim1);
		Cww_im1 = wim1_mag*wim1_mag;
		Cwu_im1 = this->dotp(wim1,ui);
		// Cuu_i = this->dotp(ui,ui);
		Cuu_i = ui_mag*ui_mag;
		Cuv_i = this->dotp(ui,vi);
	}

	// involving only i	
	// Cvv_i = this->dotp(vi,vi);
	Cvv_i = vi_mag*vi_mag;
	// Cww_i = this->dotp(wi,wi);
	Cww_i = wi_mag*wi_mag;

	// involving i+1
	if (i<N-1){
		Cwu_i = this->dotp(wi,uip1);
		// Cuu_ip1 = this->dotp(uip1,uip1);
		Cuu_ip1 = uip1_mag*uip1_mag;
		// Cvv_ip1 = this->dotp(vip1,vip1);
		Cvv_ip1 = vip1_mag*vip1_mag;
		Cuv_ip1 = this->dotp(uip1,vip1);
	}

	// calculate composite terms
	// involving i-1
	if (i>0){
		Dwu_im1 = Cww_im1*Cuu_i - Cwu_im1*Cwu_im1;
		Duv_i = Cuu_i*Cvv_i - Cuv_i*Cuv_i;
	}

	// involving i+1
	if (i<N-1){
		Dwu_i = Cww_i*Cuu_ip1 - Cwu_i*Cwu_i;
		Duv_ip1 = Cuu_ip1*Cvv_ip1 - Cuv_ip1*Cuv_ip1;
	}
	
	// get constants for force terms (K is stiffness, R is ratio of C terms)
	if (i>0){
		K1 = (eta[im1]-eta0[im1])/(sqrt(Dwu_im1));
		R1 = Cwu_im1/Cuu_i;
		K2 = (theta[i]-theta0[i])/(sqrt(Duv_i));
		R2 = (Cuv_i/Cuu_i);
	}		
	if (i<N-1){
		K3 = (eta[i]-eta0[i])/(sqrt(Dwu_i));
		R3 = (Cwu_i/Cuu_ip1);
		K4 = (theta[ip1]-theta0[ip1])/(sqrt(Duv_ip1));
		R4 = Cuv_ip1/Cuu_ip1;
	}

	// calculate force
	for (d=0; d<NDIM; d++){
		// force on N and C atoms
		fN[d] = 0;
		fC[d] = 0;

		if (i>0){
			fN[d] += kba*K1*(R1*ui[d]-wim1[d]);
			fN[d] += kba*K2*(R2*ui[d]-vi[d]);
		}
		
		if (i<N-1){
			fC[d] += kba*K3*(wi[d]-R3*uip1[d]);
			fC[d] += kba*K4*(vip1[d]-R4*uip1[d]);
		}

		// add to net force
		F[i][d] += fN[d]+fC[d];
	}	

	// get branches to N and C atoms (used below when calculating torques)
	qNx = xW[i][nid][0];
	qNy = xW[i][nid][1];
	qNz = xW[i][nid][2];

	qCx = xW[i][cid][0];
	qCy = xW[i][cid][1];
	qCz = xW[i][cid][2];


	// calculate torques due to theta angle
	TqW[i][0] += qNy * fN[2] - qNz * fN[1];
	TqW[i][1] += -qNx * fN[2] + qNz * fN[0];
	TqW[i][2] += qNx * fN[1] - qNy * fN[0];

	// calculate torques due to eta angle
	TqW[i][0] += qCy * fC[2] - qCz * fC[1];
	TqW[i][1] += -qCx * fC[2] + qCz * fC[0];
	TqW[i][2] += qCx * fC[1] - qCy * fC[0];


	// theta angle contribution
	Utheta = 0;
	if (i>0)
		Utheta = 0.5*kba*pow(theta[i]-theta0[i],2);

	// eta angle contribution
	Ueta = 0;
	if (i<N-1)
		Ueta = 0.5*kba*pow(eta[i]-eta0[i],2);

	// total potential energy
	Ubb = Utheta + Ueta;
	return Ubb;
}


// dihedral angle' force
double backbone::da_force(int i){
	// local variables
	double dx,dy,dz,Ubb;

	return Ubb;
}


// calculate dot product for term
double backbone::dotp(double a[], double b[]){
	int d = 0;
	double val = 0;

	// calculate dot product
	for (d=0; d<NDIM; d++)
		val += a[d]*b[d];

	return val;
}





/****************************
*							*
*							*
*		  PRINTING			*
*							*
*							*
*****************************/


void backbone::print_angles(){
	int i,w;

	// header
	w = 10;

	cout << "\n**\n  printing angles" << endl;
	cout << setw(w) << "theta";
	cout << setw(w) << "eta";
	cout << setw(w) << "phi_da";
	cout << setw(w) << "psi_da";
	cout << setw(w) << "omega_da";
	cout << endl;
	for (i=0; i<5*w; i++)
		cout << "=";
	cout << endl;

	// i=0
	cout << setw(w) << "/";
	cout << setw(w) << eta[0];
	cout << setw(w) << "/";
	cout << setw(w) << psi_da[0];
	cout << setw(w) << "/";
	cout << endl;
	for (i=1; i<N-1; i++){
		cout << setw(w) << theta[i];
		cout << setw(w) << eta[i];
		cout << setw(w) << phi_da[i];
		cout << setw(w) << psi_da[i];
		cout << setw(w) << omega_da[i];
		cout << endl;
	}
	cout << setw(w) << theta[N-1];
	cout << setw(w) << "/";
	cout << setw(w) << phi_da[N-1];
	cout << setw(w) << "/";
	cout << setw(w) << omega_da[N-1];
	cout << endl;	
	cout << endl;

	cout << "\n**\n leaving angle printing" << endl;

}



