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

	// initialize potential energy terms
	Ubb = 0;
	Ubl = 0;
	Uba = 0;
	Uda = 0;
	Usteric = 0;

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
	ctheta = new double [N];
	ceta = new double [N];
	cphi = new double [N];
	cpsi = new double [N];
	comega = new double [N];

	// set initial values to -1
	for (i=0; i<N; i++){
		ctheta[i] = -1.0;
		ceta[i] = -1.0;
		cphi[i] = -1.0;
		cpsi[i] = -1.0;
		comega[i] = -1.0;
	}

	// rest angles
	l0 = new double [N];
	ctheta0 = new double [N];
	ceta0 = new double [N];
	cphi0 = new double [N];
	cpsi0 = new double [N];
	comega0 = new double [N];

	// set rest lengths/angle
	cout << "initializing values of angle arrays to std values; ";	
	len0 = 1.5;					// fraction of contact distance
	this->bl0_init(len0);
	th0 = 0.6*PI;				// theta0 (everything)
	this->ba0_init(cos(th0));	
	da0 = 0;				// da0 (everything)
	this->da0_init(cos(da0));
	cout << "bl0 = " << len0 << "; ";
	cout << "ba0 = " << th0 << "; ";
	cout << "da0 = " << da0 << endl;

	cout << endl;
	cout << "test backbone constructor complete." << endl;
	cout << endl << endl;
}

// DESTRUCTOR
backbone::~backbone(){
	// free memory for values
	delete [] ctheta;
	delete [] ceta;
	delete [] cphi;
	delete [] cpsi;
	delete [] comega;

	// free memory for rest values
	delete [] ctheta0;
	delete [] ceta0;
	delete [] cphi0;
	delete [] cpsi0;
	delete [] comega0;
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
		ctheta0[i] = val;
		ceta0[i] = val;
	}
}

void backbone::da0_init(double val){
	int i;

	// loop over indices
	for (i=0; i<N; i++){
		cphi0[i] = val;
		cpsi0[i] = val;
		comega[i] = val;
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
	ceta[0] = this->get_ceta(0);				// eta depends on i+1, so has i=0 value
	cpsi[0] = this->get_cpsi(0);				// psi depends on i+1, so has i=0 value
	for (i=1; i<N-1; i++){
		ceta[i] = this->get_ceta(i);
		ctheta[i] = this->get_ctheta(i);		
		cphi[i] = this->get_cphi(i);
		cpsi[i] = this->get_cpsi(i);
		comega[i] = this->get_comega(i);
	}
	ctheta[N-1] = this->get_ctheta(N-1);			// theta depends on i-1, so has i=N-1 value
	cphi[N-1] = this->get_cphi(N-1);			// phi depends on i-1, so has i=N-1 value
	comega[N-1] = this->get_comega(N-1);		// omega depends on i-1, so has i=N-1 value
}



// GETTERS

// get theta angle for residue r
double backbone::get_ctheta(int r){
	// connection vectors
	double v1[NDIM];
	double v2[NDIM];
	double num,v1norm,v2norm,denom,val;
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
	val = -1*num/denom;

	// return value
	return val;
}

// get eta angle for residue r
double backbone::get_ceta(int r){
	// connection vector
	double v1[NDIM];
	double v2[NDIM];
	double num,v1norm,v2norm,denom,val;
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
	val = -1*num/denom;

	// return value
	return val;
}

// get phi dihedral angle for residue r (WORK IN PROGRESS)
double backbone::get_cphi(int r){
	return cos(0);
}

// get psi dihedral angle for residue r (WORK IN PROGRESS)
double backbone::get_cpsi(int r){
	return cos(0);
}

// get omega dihedral angle for residue r (WORK IN PROGRESS)
double backbone::get_comega(int r){
	return cos(0);
}



/****************************
*							*
*							*
*		 SIMULATION		    *
*							*
*							*
*****************************/



// BACKBONE: RELAX TOPOLOGY
void backbone::top_relax(int kmax){
	// local variables
	int k = 0;
	double Utol = 1e-16;
	double T0 = 1e-2;

	// energies
	Usteric = 0;
	Ubb = 0;

	// initialize velocities
	this->rand_vel_init(T0);

	// get backbone forces
	this->set_angles();
	this->print_angles();
	Ubb = this->bb_force_update(0);
	U = Ubb;

	// loop while backbone has a high potential energy
	while (k < kmax && U > Utol){
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
			this->bb_md_monitor();
			this->print_angles();
			cout << endl;
			cout << endl;
		}

		// update iterate
		k++;
	}

	// report on success or failure of loop
	if (k == kmax){
		cout << "ERROR: Loop did not relax the backbone potential energy, ending..." << endl;
		throw;
	}
	else
		cout << "Loop completed! topology of backbone sufficiently relaxed..." << endl;

	cout << "Loop completed! topology of backbone sufficiently relaxed..." << endl;
	/*
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
			this->bb_md_monitor();
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
		*/
}

void backbone::bb_free_md(double tmp0, int NT, int nnu) {
	int t;

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
			this->bb_md_monitor();
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

	// reset Usterictmp
	double Usterictmp = 0;
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
								Usterictmp += (ep / 2) * pow(1 - da / rij, 2);
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
								Usterictmp += (ep / 2) * pow(1 - da / rij, 2);
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
	return Usterictmp;
}


// BACKBONE FORCE UPDATE

// choice can choose which interaction to turn on:
// choice = 0 -> only BL force
// choice = 1 -> only BA force
// choice = 2 -> only DA force
double backbone::bb_force_update(int choice){
	// local variables
	int i,ip1;
	double Ubbtmp = 0;
	double Utmp = 0;

	// bl, ba, da energies
	Ubl = 0;
	Uba = 0;
	Uda = 0;

	// loop over particles, calculate forces due to backbone connections
	for (i=0; i<N; i++){
		// bond length force (done pairwise)
		if (i > 0){
			Utmp = this->bl_force(i);
			Ubl += Utmp;
			Ubbtmp += Utmp;
		}

		// bond angle force
		if (choice > 0){
			if (i < N-1){
				Utmp = this->ba_force(i);
				Uba += Utmp;
				Ubbtmp += Utmp;
			}			
		}

		// // dihedral angle force
		if (choice > 1){
			Utmp = this->bl_force(i);
			Uda += Utmp;
			Ubbtmp += Utmp;
		}
	}
	return Ubbtmp;
}

// bond length force
double backbone::bl_force(int i){
	// local variables
	double dij[NDIM];
	double vfij[NDIM];
	double rij,fij,r0,Ubbtmp;
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
	Ubbtmp = 0.5*kbl*pow(rij-r0,2);

	return Ubbtmp;
}


// bond angle force

// forces determined by forces on Ci
double backbone::ba_force(int i){
	// indx variables
	int d;
	int ip1 = i+1;

	// force and potential variables
	double f[NDIM];							// f vector, force on Ci atom due to bond-angle force
	double Utheta, Ueta, Ubbtmp;			// potential energies

	// constants in force term
	double delceta, delctheta;				// distances to rest angle
	double Cww_i, Cuu_ip1, Cvv_ip1;			// dot product terms
	double Mwu_i, Muv_ip1;					// sqrt terms
	double K1, K2, K3;						// constants in force term (see latex doc, section 3)
	double wi_mag, uip1_mag, vip1_mag;		// vector magnitudes
	double wi[NDIM];						// wi vector, Cai -> Ci
	double uip1[NDIM];						// uip1 vector, Ci -> Nip1
	double vip1[NDIM];						// vip1 vector, Nip1 -> Caip1	

	// branch vectors for torque calc
	double qix, qiy, qiz;
	double qip1x, qip1y, qip1z;

	// calculate distances from rest angle
	delceta = ceta[i] - ceta0[i];
	delctheta = ctheta[ip1] - ctheta0[ip1];
	
	// calculate vectors and magnitudes
	wi_mag = this->get_atomic_distance(i,i,caid,cid,wi);
	uip1_mag = this->get_atomic_distance(i,ip1,cid,nid,uip1);
	vip1_mag = this->get_atomic_distance(ip1,ip1,nid,caid,vip1);

	// calculate dot product terms
	Cww_i = wi_mag*wi_mag;
	Cuu_ip1 = uip1_mag*uip1_mag;
	Cvv_ip1 = vip1_mag*vip1_mag;

	// calculate sqrt terms
	Mwu_i = 1/sqrt(Cww_i*Cuu_ip1);
	Muv_ip1 = 1/sqrt(Cuu_ip1*Cvv_ip1);

	// calculate larger constants
	K1 = (delceta*ceta[i] + delctheta*ctheta[ip1])/Cuu_ip1;
	K2 = delceta*Mwu_i;
	K3 = delctheta*Muv_ip1;

	// calculate force
	for (d=0; d<NDIM; d++){
		// force on Ci (and therefore Ni+1)
		f[d] = kba*(K1*uip1[d] + K2*wi[d] + K3*vip1[d]);	

		// add to net force on i (due to atom Ci)
		F[i][d] += f[d];

		// by reciprocity, add to net force on ip1 (due to atom Nip1)
		F[ip1][d] -= f[d];
	}	

	// get branches to Ci and Nip1 atoms (used below when calculating torques)
	qix = xW[i][cid][0];
	qiy = xW[i][cid][1];
	qiz = xW[i][cid][2];

	qip1x = xW[ip1][nid][0];
	qip1y = xW[ip1][nid][1];
	qip1z = xW[ip1][nid][2];


	// calculate torques due to force on Ci
	TqW[i][0] += qiy * f[2] - qiz * f[1];
	TqW[i][1] += -qix * f[2] + qiz * f[0];
	TqW[i][2] += qix * f[1] - qiy * f[0];

	// calculate torques due to force on Nip1
	TqW[ip1][0] -= qip1y * f[2] - qip1z * f[1];
	TqW[ip1][1] -= -qip1x * f[2] + qip1z * f[0];
	TqW[ip1][2] -= qip1x * f[1] - qip1y * f[0];

	// eta angle contribution
	Ueta = 0.5*kba*delceta*delceta;

	// theta angle contribution
	Utheta = 0.5*kba*delctheta*delctheta;	

	// total potential energy
	Ubbtmp = Utheta + Ueta;
	return Ubbtmp;
}


// dihedral angle' force
double backbone::da_force(int i){
	// local variables
	double dx,dy,dz,Ubbtmp;

	return Ubbtmp;
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
	cout << setw(w) << "ctheta";
	cout << setw(w) << "ceta";
	cout << setw(w) << "cphi";
	cout << setw(w) << "cpsi";
	cout << setw(w) << "comega";
	cout << endl;
	for (i=0; i<5*w; i++)
		cout << "=";
	cout << endl;

	// i=0
	cout << setw(w) << "/";
	cout << setw(w) << ceta[0];
	cout << setw(w) << "/";
	cout << setw(w) << cpsi[0];
	cout << setw(w) << "/";
	cout << endl;
	for (i=1; i<N-1; i++){
		cout << setw(w) << ctheta[i];
		cout << setw(w) << ceta[i];
		cout << setw(w) << cphi[i];
		cout << setw(w) << cpsi[i];
		cout << setw(w) << comega[i];
		cout << endl;
	}
	cout << setw(w) << ctheta[N-1];
	cout << setw(w) << "/";
	cout << setw(w) << cphi[N-1];
	cout << setw(w) << "/";
	cout << setw(w) << comega[N-1];
	cout << endl;	
	cout << endl;

	cout << "\n**\n  printing rest cosine angles" << endl;
	cout << setw(w) << "ctheta0";
	cout << setw(w) << "ceta0";
	cout << setw(w) << "cphi0";
	cout << setw(w) << "cpsi0";
	cout << setw(w) << "comega0";
	cout << endl;
	for (i=0; i<5*w; i++)
		cout << "=";
	cout << endl;

	// i=0
	cout << setw(w) << "/";
	cout << setw(w) << ceta0[0];
	cout << setw(w) << "/";
	cout << setw(w) << cpsi0[0];
	cout << setw(w) << "/";
	cout << endl;
	for (i=1; i<N-1; i++){
		cout << setw(w) << ctheta0[i];
		cout << setw(w) << ceta0[i];
		cout << setw(w) << cphi0[i];
		cout << setw(w) << cpsi0[i];
		cout << setw(w) << comega0[i];
		cout << endl;
	}
	cout << setw(w) << ctheta0[N-1];
	cout << setw(w) << "/";
	cout << setw(w) << cphi0[N-1];
	cout << setw(w) << "/";
	cout << setw(w) << comega0[N-1];
	cout << endl;	
	cout << endl;

	cout << "\n**\n  printing angle distance to rest" << endl;
	cout << setw(w) << "ctheta";
	cout << setw(w) << "ceta";
	cout << setw(w) << "cphi";
	cout << setw(w) << "cpsi";
	cout << setw(w) << "comega";
	cout << endl;
	for (i=0; i<5*w; i++)
		cout << "=";
	cout << endl;

	// i=0
	cout << setw(w) << "/";
	cout << setw(w) << ceta[0]-ceta0[0];
	cout << setw(w) << "/";
	cout << setw(w) << cpsi[0]-cpsi0[0];
	cout << setw(w) << "/";
	cout << endl;
	for (i=1; i<N-1; i++){
		cout << setw(w) << ctheta[i]-ctheta0[i];
		cout << setw(w) << ceta[i]-ceta0[i];
		cout << setw(w) << cphi[i]-cphi0[i];
		cout << setw(w) << cpsi[i]-cpsi0[i];
		cout << setw(w) << comega[i]-comega0[i];
		cout << endl;
	}
	cout << setw(w) << ctheta[N-1]-ctheta0[N-1];
	cout << setw(w) << "/";
	cout << setw(w) << cphi[N-1]-cphi0[N-1];
	cout << setw(w) << "/";
	cout << setw(w) << comega[N-1]-comega0[N-1];
	cout << endl;	
	cout << endl;

	cout << "\n**\n leaving angle printing" << endl;
}


void backbone::bb_md_monitor() {
	cout << "** Packing:" << endl;
	cout << "N = " << N << endl;
	cout << "phi = " << phi << endl;
	cout << endl;
	cout << "** Energy:" << endl;
	cout << "U = " << U/N << endl;
	cout << "K = " << K/N << endl;
	cout << "Krot = " << Krot/N << endl;
	cout << "Ktrans = " << K/N - Krot/N << endl;
	cout << "E = " << U/N + K/N << endl;
	cout << endl;
	cout << "** Contacts:" << endl;
	cout << "sum c = " << this->get_c_sum() << endl;
	cout << "sum ac = " << 0.5 * (this->get_ac_sum()) << endl;
	cout << "niso max = " << DOF*N - NDIM + 1 << endl;
	cout << endl;
	cout << "** FIRE:" << endl;
	cout << "alpha = " << alpha << endl;
	cout << "dt = " << dt << endl;
	cout << "dtmax = " << dtmax << endl;
	cout << endl;

	// output to energy file if open
	if (enobj.is_open()) {
		cout << "Printing ENERGY" << endl;
		enobj << setw(12) << U/N;
		enobj << setw(12) << Ubb/N;
		enobj << setw(12) << Ubl/N;
		enobj << setw(12) << Uba/N;
		enobj << setw(12) << Uda/N;
		enobj << setw(12) << Usteric/N;
		enobj << setw(12) << K/N;
		enobj << setw(12) << Krot/N;
		enobj << setw(12) << U/N + K/N;
		enobj << endl;
	}	

	// output to xyz file if open
	if (xyzobj.is_open()) {
		cout << "Printing XYZ" << endl;
		this->rigidbody_xyz();
		cout << endl;
	}
}

