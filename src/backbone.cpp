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
	caid = 0;
	cid = 0;

	// set spring constants to 1
	kbl = 1.0;
	kba = 1.0;
	kda = 1.0;

	// allocate memory
	cout << "allocating memory to angle and rest angle arrays..." << endl;

	// arrays of values
	theta = new double [N-1];
	eta = new double [N-1];
	phi_da = new double [N-1];
	psi_da = new double [N-1];
	omega_da = new double [N-1];

	for (i=0; i<N-1; i++){
		theta[i] = 0;
		eta[i] = 0;
		phi_da[i] = 0;
		psi_da[i] = 0;
		omega_da[i] = 0;
	}

	// rest angles
	l0 = new double [N-1];
	theta0 = new double [N-1];
	eta0 = new double [N-1];
	phi0_da = new double [N-1];
	psi0_da = new double [N-1];
	omega0_da = new double [N-1];

	// set rest lengths/angle
	cout << "initializing values of angle arrays to std values; ";
	cout << "bl0 = " << len0 << "; ";
	cout << "ba0 = " << th0 << "; ";
	cout << "da0 = " << da0 << endl;
	len0 = 0.8;					// fraction of contact distance
	this->bl0_init(len0);
	th0 = 0.2*PI;				// theta0 (everything)
	this->ba0_init(th0);	
	da0 = 0.2*PI;				// da0 (everything)
	this->da0_init(da0);

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
	for (i=0; i<N-1; i++)
		l0[i] = val;
}

void backbone::ba0_init(double val){
	int i;

	// loop over indices
	for (i=0; i<N-1; i++){
		theta0[i] = val;
		eta0[i] = val;
	}
}

void backbone::da0_init(double val){
	int i;

	// loop over indices
	for (i=0; i<N-1; i++){
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
	v1norm = 0;
	v2norm = 0;
	for (d=0; d<NDIM; d++){
		// define vector components
		v1[d] = xW[r][nid][d]-xW[rm1][cid][d];
		v2[d] = xW[r][caid][d]-xW[r][nid][d];

		// define norms
		v1norm += v1[d]*v1[d];
		v2norm += v2[d]*v2[d];
	}

	// numerator (use definition of C)
	num = this->C_prod(v1,v2);

	// denominator
	denom = sqrt(v1norm*v2norm);

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
	v1norm = 0;
	v2norm = 0;
	for (d=0; d<NDIM; d++){
		// define vector components
		v1[d] = xW[r][cid][d]-xW[r][caid][d];
		v2[d] = xW[rp1][nid][d]-xW[r][cid][d];		

		// define norms
		v1norm += v1[d]*v1[d];
		v2norm += v2[d]*v2[d];
	}

	// numerator (use definition of C)
	num = this->C_prod(v1,v2);

	// denominator
	denom = sqrt(v1norm*v2norm);

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
	int kmax = 1e3;
	double Ubb = 0;
	double Utol = 1e-16;
	double T0 = 1e-2;

	// initialize velocities
	this->rand_vel_init(T0);

	// get backbone forces
	this->set_angles();
	this->print_angles();
	Ubb = this->bb_force_update();
	k = kmax;

	// loop while backbone has a high potential energy
	while (Ubb > Utol && k < kmax){
		// advance quaternions, positions
		this->verlet_first();

		// update angles
		this->set_angles();

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




// BACKBONE FORCE UPDATE

double backbone::bb_force_update(){
	// local variables
	int i,ip1;
	double Ubb = 0;

	cout << "\n**\n  In bb force update!" << endl;

	// loop over particles, calculate forces due to backbone connections
	for (i=0; i<N; i++){
		// bond length force (done pairwise)
		if (i > 0)
			Ubb = Ubb + this->bl_force(i);
		cout << "i = " << i << ", after bl Ubb = " << Ubb << "; ";

		// bond angle force
		Ubb = Ubb + this->ba_force(i);
		cout << "after ba, Ubb = " << Ubb << endl;

		// // dihedral angle force
		// U = U + this->da_force(i);
	}

	cout << "  Leaving bb force update!\n**" << endl << endl;

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
	double ftheta[NDIM];				// force vector on particle i due to theta angle
	double feta[NDIM];					// force vector on particle i due to eta angle
	double qNx,qNy,qNz,qCx,qCy,qCz;		// branch vectors to atomic coordinates

	// connection vectors
	double wim1[NDIM];
	double ui[NDIM];
	double vi[NDIM];
	double wi[NDIM];
	double uip1[NDIM];
	double vip1[NDIM];

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
	for (d=0; d<NDIM; d++){
		// vectors with prevous atomic positions (do not exist for i = 0)
		if (i>0){
			wim1[d] = xW[im1][cid][d] - xW[im1][caid][d];
			ui[d] = xW[i][nid][d] - xW[im1][cid][d];			
		}

		// vectors with atoms of residue i
		vi[d] = xW[i][caid][d] - xW[i][nid][d];
		wi[d] = xW[i][cid][d] - xW[i][caid][d];

		// vectors with next atom positions (do not exist for i = N-1)
		if (i<N-1){
			uip1[d] = xW[ip1][nid][d] - xW[i][cid][d];
			vip1[d] = xW[ip1][caid][d] - xW[ip1][nid][d];
		}
	}

	// calculate dot product terms
	// involving i-1
	if (i>0){
		Cww_im1 = this->C_prod(wim1,wim1);
		Cwu_im1 = this->C_prod(wim1,ui);
		Cuu_i = this->C_prod(ui,ui);
		Cuv_i = this->C_prod(ui,vi);
	}

	// involving only i	
	Cvv_i = this->C_prod(vi,vi);	
	Cww_i = this->C_prod(wi,wi);

	// involving i+1
	if (i<N-1){
		Cwu_i = this->C_prod(wi,uip1);
		Cuu_ip1 = this->C_prod(uip1,uip1);
		Cvv_ip1 = this->C_prod(vip1,vip1);
		Cuv_ip1 = this->C_prod(uip1,vip1);
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
		K1 = (eta[im1]-eta0[im1])/sqrt(Dwu_im1);			
		R1 = Cwu_im1/Cuu_i;
		K2 = (theta[i]-theta0[i])/sqrt(Duv_i);
		R2 = Cuv_i/Cuu_i - 2;
	}		
	if (i<N-1){
		K3 = (eta[i]-eta0[i])/sqrt(Dwu_i);
		R3 = 2 - Cwu_i/Cuu_ip1;
		K4 = (theta[ip1]-theta0[ip1])/sqrt(Duv_ip1);
		R4 = Cuv_ip1/Cuu_ip1;
	}

	// calculate force
	for (d=0; d<NDIM; d++){
		// force due to theta angle
		ftheta[d] = 0;
		if (i>0){
			ftheta[d] += kba*K1*(R1*ui[d]-wim1[d]);
			ftheta[d] += kba*K2*(R2*ui[d]+vi[d]);
		}

		// force due to eta angle
		feta[d] = 0;
		if (i<N-1){
			feta[d] += kba*K3*(R3*ui[d]-wi[d]);
			feta[d] += kba*K4*(vip1[d]-R4*uip1[d]);
		}

		// add to net force
		F[i][d] += ftheta[d]+feta[d];
	}	

	// get branches to N and C atoms (used below when calculating torques)
	qNx = xW[i][nid][0];
	qNy = xW[i][nid][1];
	qNz = xW[i][nid][2];

	qCx = xW[i][cid][0];
	qCy = xW[i][cid][1];
	qCz = xW[i][cid][2];


	// calculate torques due to theta angle
	TqW[i][0] += qNy * ftheta[2] - qNz * ftheta[1];
	TqW[i][1] += -qNx * ftheta[2] + qNz * ftheta[0];
	TqW[i][2] += qNx * ftheta[1] - qNy * ftheta[0];

	// calculate torques due to eta angle
	TqW[i][0] += qCy * feta[2] - qCz * feta[1];
	TqW[i][1] += -qCx * feta[2] + qCz * feta[0];
	TqW[i][2] += qCx * feta[1] - qCy * feta[0];


	// theta angle contribution
	Utheta = 0;
	if (i>0)
		Utheta = 0.5*kba*(theta[i]-theta0[i]);

	// eta angle contribution
	Ueta = 0;
	if (i<N-1)
		Ueta = 0.5*kba*(eta[i]-eta0[i]);

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





// calculate C term
double backbone::C_prod(double a[], double b[]){
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

	for (i=0; i<N-1; i++){
		cout << setw(w) << theta[i];
		cout << setw(w) << eta[i];
		cout << setw(w) << phi_da[i];
		cout << setw(w) << psi_da[i];
		cout << setw(w) << omega_da[i];
		cout << endl;
	}	
	cout << endl;

	cout << "\n**\n leaving angle printing" << endl;

}



