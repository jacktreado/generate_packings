/*

	Methods implementation 
	for packing class

	BY Jack Treado

*/

#include "packing.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

const double PI = 3.1415926;




/* 
==================================

		  INITIALIZATION		 

================================== 
*/



void packing::mxwl_vel_init (double t)
/* This subroutine assigns velocities according to the Maxwell
 distribution. N is the total number of velocity components and
 M is the total number of degrees of freedom. T is the system
 temperature in the reduced units.  Copyright (c) Tao Pang 1997. */
{
	int i,d;
    double v1;
    double ek,vs;

    void grnf(double&);
    
    /* seed random number generator */
    srandom((unsigned int) seed);

    /* Assign a Gaussian distribution to each velocity component */   
    for (i=0; i<N; i++){
    	for (d=0; d<NDIM; d++){
    		grnf(v1);
        	v[i][d] = v1;
    	}        
    }
    
    /* Scale the velocity to satisfy the equipartition theorem */
    ek = 0;
    for (i=0; i<N; i++){
    	for (d=0; d<NDIM; d++)
    		ek += 0.5*m[i]*v[i][d]*v[i][d];
    }
    ek /= N;

    vs = sqrt(ek/t);
    for (i=0; i<N; i++){
    	for (d=0; d<NDIM; d++)
    		v[i][d] /= vs;
    }
}

void grnf (double& x)
/* Two Gaussian random numbers generated from two uniform
 random numbers.  Copyright (c) Tao Pang 1997. */
{
    double pi,r1,r2,r;
    
    pi =  4*atan(1.0);
    r = (double)random()/(RAND_MAX+1.0);
    r1 = -log(1-r);
    r = (double)random()/(RAND_MAX+1.0);
    r2 =  2*pi*r;
    r1 =  sqrt(2*r1);
    x  = r1*cos(r2);
}

void packing::rand_vel_init(double tmp0){
	int i,d;
	double ek,vs;
	double* pmean;
	srand48(seed);

	pmean = new double[NDIM];
	for (d=0; d<NDIM; d++)
		pmean[d] = 0.0;

	for (i=0; i<N; i++){
    	for (d=0; d<NDIM; d++){
        	v[i][d] = drand48();
        	pmean[d] += m[i]*v[i][d];        	
    	}
    }

    for (d=0; d<NDIM; d++)
		pmean[d] /= N;

	for (i=0; i<N; i++){
    	for (d=0; d<NDIM; d++){
        	v[i][d] -= pmean[d]/m[i];
        	ek = 0.0;
        	ek += 0.5*m[i]*v[i][d]*v[i][d];
    	}
    }

    vs = sqrt(ek/tmp0);
    for (i=0; i<N; i++){
    	for (d=0; d<NDIM; d++)
        	v[i][d] /= vs;
    }

    delete [] pmean;
}






/* 
==================================

		Molecular Dynamics		 

================================== 
*/

// harmonic spring interactions
void packing::md_evol_hs(double tmp0, int NT){

	// initialize particle velocities
	cout << "* initializing MD velocities..." << endl;
	this->rand_vel_init(tmp0);

	// loop over time
	int t;
	for (t=0; t<NT; t++){
		// first md step
		this->pos_update();

		// force update
		this->hs_force();

		// vel update
		this->vel_update();

		// print info
		if (t % plotskip == 0){
			cout << "~~~~ only MD ~~~~" << endl;
			this->md_monitor(t,0,-1,-1);
			if (xyzobj.is_open())
				this->print_xyz();
			cout << endl << endl;
		}
	}

	cout << "MD completed!" << endl;
}

// update particle positions using velocity-verlet
void packing::pos_update(){
	int i,d;
	double pos;

	// update kinetic energy
	K = 0;

	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			pos = x[i][d];

			// update positions
			pos += dt*v[i][d] + 0.5*dt*dt*aold[i][d];

			// check periodic boundaries
			if (pos < 0)
				pos += L[d];
			else if (pos > L[d])
				pos -= L[d];				

			x[i][d] = pos;

			// update K
			K += 0.5*m[i]*v[i][d]*v[i][d];

			// reset forces
			F[i][d] = 0;
		}
		// reset contacts
		pc[i] = 0;			
	}
}


/*	FORCE UPDATES  */
// HS: harmonic spring
// HZS: hertzian spring
// RLJ: repulsive lennard jones
// ES: electrostatics

double packing::hs(double sij, double xij){
	double f;

	// calc pairwise force, update potential
	f = -(ep/(sij*sij))*(sij/xij-1);

	return f;
}

void packing::hs_force(){
	int i,j,d,cind;
	double sij, dx;
	double xij[NDIM];

	// reset U
	U = 0;

	for (i=0; i<N; i++){
		for (j=i+1; j<N; j++){
			// contact matrix index
			cind = N*i + j - ((i+1)*(i+2))/2;

			// get contact distance sij
			sij = r[j] + r[i];

			// get distance between particles
			dx = this->get_distance(i,j,xij);

			// particles in contact, calculate force
			if (dx < sij){				
				for (d=0; d<NDIM; d++){
					F[i][d] += (this->hs(sij,dx))*xij[d];
					F[j][d] -= (this->hs(sij,dx))*xij[d];
				}

				// update contact matrix			
				c[cind] = 1;

				// update contact list
				pc[i]++;
				pc[j]++;				

				// update potential energy
				U += (ep/2)*pow(1-dx/sij,2);
			}
			// else, no force added and no contact			
			else
				c[cind] = 0;
		}
	}
}

void packing::hs_force_nn(){
	int i,j,jj,NN,cind,d;
	double dx,sij; 
	double xij[NDIM];

	// reset U
	U = 0;

	for (i=0; i<N; i++){
		NN = neighborlist[i].size();
		for (jj=0; jj<NN; jj++){
			// get neighbor from neighborlist
			j = neighborlist[i].at(jj);

			// contact matrix index
			cind = N*i + j - ((i+1)*(i+2))/2;

			// get contact distance sij
			sij = r[j] + r[i];

			// get distance between particles
			dx = this->get_distance(i,j,xij);

			// particles in contact, calculate force
			if (dx < sij){				

				for (d=0; d<NDIM; d++){
					F[i][d] += (this->hs(sij,dx))*xij[d];
					// F[j][d] -= (this->hs(sij,dx))*xij[d];
				}

				// update contact matrix			
				c[cind] = 1;

				// update contact list
				pc[i]++;
				// pc[j]++;				

				// update potential energy
				U += (ep/2)*pow(1-dx/sij,2);
			}
			// else, no force added and no contact			
			else
				c[cind] = 0;
		}
	}
}


void packing::vel_update(){
	int i,d;
	double vel, anew;
	double* vmean;

	// vmean = new double[NDIM];
	// for (d=0; d<NDIM; d++)
	// 	vmean[d] = 0.0;

	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			vel = v[i][d];

			// update velocity
			anew = F[i][d]/m[i];
			vel += 0.5*dt*(anew + aold[i][d]);
			// vmean[d] += vel;
			v[i][d] = vel;
			aold[i][d] = anew;
		}
	}

	// for (d=0; d<NDIM; d++)
	// 	vmean[d] /= N;

	// for (i=0; i<N; i++){
	// 	for (d=0; d<NDIM; d++)
	// 		v[i][d] -= vmean[d];
	// }

	// delete [] vmean;
}


// remove particles with unconstranted dof
int packing::rmv_rattlers(int& krcrs) {	
	int i,j,ci,cj,r,nr,nm;	

	// monitor recursion depth
	krcrs++;

	// number of rattlers
	nr = 0;

	// number of "marginal" rattlers to be removed
	nm = 0;

	// loop over rows, eliminate contacts to rattlers
	for (i=0; i<N; i++){
		// start with r=0 contacts
		r = 0;

		// count contacts
		for (j=0; j<N; j++){			
			if (j < i){
				ci = N*j + i - ((j+1)*(j+2))/2;
				r += c[ci];
			}
			else if (j > i){
				cj = N*i + j - ((i+1)*(i+2))/2;	// mapping from matrix space to sub matrix space	
				r += c[cj];
			}		
			else
				continue;
		}

		// remove from network if r <= DOF, delete contacts
		if (r <= DOF){
			// increment # of rattlers
			nr++;

			// alter pc
			pc[i] = 0;

			// if in contact, remove contacts
			if (r > 0){	
				nm++;				
				for (j=0; j<N; j++){			
					if (j < i){
						ci = N*j + i - ((j+1)*(j+2))/2;
						if (c[ci] > 0){
							c[ci] = 0;
							pc[j]--;
						}						
					}
					else if (j > i){
						cj = N*i + j - ((i+1)*(i+2))/2;	// mapping from matrix space to sub matrix space	
						if (c[cj] > 0){
							c[cj] = 0;
							pc[j]--;
						}
					}	
					else
						continue;
				}
			}
		}
	}

	if (krcrs > 100){
		cout << "max recursive depth reached, be wary of output" << endl;
		return -1;
	}

	// recursively check once rattlers are removed
	if (nm == 0)
		return nr;
	else
		return rmv_rattlers(krcrs);
}



/* 
==================================

		 Time Evolution		 

================================== 
*/

// generate jammed packing using fire
void packing::jamming_finder(double tend, double dphi, double Utol, double Ktol) {
	// local variables
	int t,kr;
	double dphi0,tmp0,phiL,phiH;
	int check_rattlers;

	// set initial temperature
	tmp0 = 1.0;

	// initialize particle velocities
	cout << "* initializing velocities in jamming_finder()..." << endl;
	this->mxwl_vel_init(tmp0);

	// initialize variables
	kr = 0;
	dphi0 = dphi;
	phiL = -1;
	phiH = -1;
	check_rattlers = 0;

	// constant energy checking
	double Uold,dU,dUtol; 
	int epc,epcN,epconst;
	epconst = 0;
	Uold = 0;
	dU = 0;
	dUtol = 1e-6;
	epc = 0;
	epcN = 5e3;

	// get number of time steps
	int NT;
	NT = round(tend/dt);
	dt = tend/NT;

	cout << "===== BEGINNING TIME LOOP ======" << endl;
	cout << " NT = " << NT << endl;
	cout << " dt = " << dt << endl;
	cout << " tend = " << tend << endl;
	cout << "================================" << endl << endl;

	// loop over time
	for (t=0; t<NT; t++){		
		// first md time step
		this->pos_update();

		// force update
		this->hs_force();

		// vel update
		this->vel_update();

		// check for nan
		// this->catch_nans(t);

		// fire update (inertial relaxer)
		this->fire();

		// check for rattlers if root search is on
		if (check_rattlers){
			kr = 0;
			nr = this->rmv_rattlers(kr);
		}		
		else{
			nr = 0;
		}

		if (t % plotskip == 0){
			this->md_monitor(t,nr,phiH,phiL);
			cout << endl;
			cout << "** Utol = " << Utol << endl;
			cout << "** Ktol = " << Ktol << endl;
			cout << "** check_rattlers = " << check_rattlers << endl;
			cout << endl << endl;
			if (xyzobj.is_open())
				this->print_xyz();
		}

		if (phi > 1.0){
			this->print_vars();
			throw "phi > 1, ending...\n";
		}

		// check for constant potential energy
		dU = abs(Uold-U);
		if (dU < dUtol){
			epc++;
			if (epc > epcN)
				epconst = 1;
		}
		else{
			epconst = 0;
			epc = 0;
		}
		Uold = U;

		if (epconst == 1)
			check_rattlers = 1;

		this->root_search(phiH,phiL,check_rattlers,epconst,nr,dphi0,Ktol,Utol,t);

		if (isjammed == 1)
			break;
	}
}


// generate jammed packing using fire
void packing::jamming_finder_nn(double tend, double dphi, double Utol, double Ktol) {
	// local variables
	int t,kr,csum,niso,nbb;
	double dphi0,tmp0,phiL,phiH;
	bool oc,uc,gr,jammed,marginal;
	int check_rattlers;

	// set initial temperature
	tmp0 = 5.0;

	// initialize particle velocities
	cout << "* initializing velocities in jamming_finder()..." << endl;
	this->mxwl_vel_init(tmp0);

	// initialize variables
	kr = 0;
	dphi0 = dphi;
	phiL = -1;
	phiH = -1;
	check_rattlers = 0;

	// constant energy checking
	double Uold,dU,dUtol; 
	int epc,epcN,epconst;
	epconst = 0;
	Uold = 0;
	dU = 0;
	dUtol = 1e-6;
	epc = 0;
	epcN = 1e3;

	// get number of time steps
	int NT;
	NT = round(tend/dt);
	dt = tend/NT;

	// get cell neighbor list and positions
	this->get_cell_neighbors();
	this->get_cell_positions();

	cout << "===== BEGINNING TIME LOOP ======" << endl;
	cout << "NT = " << NT << endl;
	cout << "dt = " << dt << endl;
	cout << "tend = " << tend << endl;

	// loop over time
	for (t=0; t<NT; t++){	
		// update nearest neighbor lists
		if (t % nnupdate == 0){
			// cout << "updating cell list" << endl;
			this->update_cell();
			// cout << "updating neighborlist" << endl;
			this->update_neighborlist();
		}

		if (t % plotskip == 0){
			// cout << "printing??" << endl;
			this->md_monitor(t,nr,phiH,phiL);
			cout << endl;
			cout << "** check_rattlers = " << check_rattlers << endl;
			cout << "** epconst = " << epconst << endl;
			if (xyzobj.is_open()){
				// cout << "PRINTING NL INFO..." << endl;
				// this->print_cell();
				// this->print_neighborlist();
				this->print_nl_xyz();									
			}
			cout << endl << endl;
		}	

		// first md time step
		this->pos_update();

		// force update
		this->hs_force_nn();

		// vel update
		this->vel_update();

		// check for nan
		this->catch_nans(t);

		// fire update (inertial relaxer)
		if (t > 500)
			this->fire();

		// check for rattlers if root search is on
		if (check_rattlers){
			kr = 0;
			nr = this->rmv_rattlers(kr);
		}		
		else
			nr = 0;

		nbb = N - nr;
		niso = DOF*nbb - NDIM + 1;
		csum = this->get_c_sum();		

		if (phi > 1.0){
			this->print_vars();
			throw "phi > 1, ending...\n";
		}

		// check for constant potential energy
		dU = abs(Uold-U);
		if (dU < dUtol){
			epc++;
			if (epc > epcN)
				epconst = 1;
		}
		else{
			epconst = 0;
			epc = 0;
		}
		Uold = U;

		
		if (t > 500)
			this->root_search(phiH,phiL,check_rattlers,epconst,nr,dphi0,Ktol,Utol,t);
		if (isjammed == 1)
			break;
	}
}


void packing::root_search(double& phiH, double& phiL, int& check_rattlers, int epconst, int nr, double dphi0, double Ktol, double Utol, int t){
	/* 
		GROWTH ALGORITHM

		1. If U < Utol grow
		2. If U > 2*Utol & K < Ktol, stop, shrink by dphi0/2
			a. If U < Utol & K < Ktol, relaxed, so grow again
			b. If U > 2*Utol & K < Ktol, still overcompressed, so shrink by dphi0/2
		3. Jammed when K < Ktol & Utol < U < 2Utol
	*/

	int nbb,niso,csum;
	bool gr,oc,uc,marginal,jammed;
	double dphi = dphi0;

	gr = 0;
	oc = 0;
	uc = 0;
	marginal = 0;
	jammed = 0;

	nbb = N - nr;
	niso = DOF*nbb - NDIM + 1;
	csum = this->get_c_sum();

	gr = (U < Utol);
	oc = (U > Utol && K < Ktol && csum >= niso && epconst == 1);
	uc = (U < Utol);
	marginal = (K < Ktol && nr == N && epconst);
	jammed = (U > Utol && U < 2*Utol && K < Ktol && csum == niso && nr < N);


	if (phiH < 0){
		if (gr){
			check_rattlers = 0;
			this->scale_sys(dphi);
		}
		else if (oc){			
			phiH = phi;
			dphi = -0.05*dphi0;
			check_rattlers = 1;

			cout << endl;
			cout << "phiH 1st set at nt = " << t << endl;
			cout << "phiH = " << phi << endl;
			cout << "phiL = " << phiL << endl;
			cout << "dphi = " << dphi << endl;
			cout << "old phi = " << phi << endl;
			this->scale_sys(dphi);
			cout << "new phi = " << phi << endl;
			cout << "entering root search..." << endl;
			cout << endl;
		}
	}
	else{
		if (phiL < 0){

			// if still overcompressed, decrease again
			if (oc){
				phiH = phi;
				dphi = -0.05*dphi0;

				cout << endl;
				cout << "still overcompressed..." << endl;
				cout << "phiH set at nt = " << t << endl;
				cout << "phiH = " << phi << endl;
				cout << "phiL = " << phiL << endl;
				cout << "dphi = " << dphi << endl;					
				cout << "old phi = " << phi << endl;
				this->scale_sys(dphi);
				cout << "new phi = " << phi << endl;					
				cout << "continuing root search..." << endl;
			}

			// if undercompressed, set phiL, root search
			if (uc){
				phiL = phi;
				dphi = 0.5*(phiH + phiL) - phi;

				cout << endl;
				cout << "relaxation found!" << endl;
				cout << "phiL set at nt = " << t << endl;
				cout << "phiH = " << phi << endl;
				cout << "phiL = " << phiL << endl;
				cout << "dphi = " << dphi << endl;					
				cout << "old phi = " << phi << endl;
				this->scale_sys(dphi);
				cout << "new phi = " << phi << endl;
				cout << "root search interval found, growing..." << endl;
			}

			// if marginal, grow, reset root search
			if (marginal){
				phiL = -1;
				phiH = -1;
				dphi = (1+0.05*drand48())*dphi0;

				cout << endl;
				cout << "marginal state found..." << endl;
				cout << "root search reset at nt = " << t << endl;
				cout << "phiH = " << phi << endl;
				cout << "phiL = " << phiL << endl;
				cout << "dphi = " << dphi << endl;
				cout << "old phi = " << phi << endl;
				this->scale_sys(dphi);
				cout << "new phi = " << phi << endl;
				cout << "marginal, so trying root search again..." << endl;
			}

			if (jammed){
				phiL = 0.9*phi;
				dphi = 0.5*(phiH + phiL) - phi;

				cout << endl;
				cout << "almost jammed found!" << endl;
				cout << "phiL set at nt = " << t << endl;
				cout << "phiH = " << phi << endl;
				cout << "phiL = " << phiL << endl;
				cout << "dphi = " << dphi << endl;					
				cout << "old phi = " << phi << endl;
				this->scale_sys(dphi);
				cout << "new phi = " << phi << endl;
				cout << "root search interval found, growing..." << endl;
			}

		}
		else{

			// if overcompressed, root search down
			if (oc){
				phiH = phi;
				dphi = 0.5*(phiH + phiL) - phi;

				cout << endl;
				cout << "overcompressed state found!" << endl;
				cout << "phiH set at nt = " << t << endl;
				cout << "phiH = " << phiH << endl;
				cout << "phiL = " << phiL << endl;
				cout << "dphi = " << dphi << endl;
				cout << "old phi = " << phi << endl;
				this->scale_sys(dphi);
				cout << "new phi = " << phi << endl;
				cout << "decreasing root search guess..." << endl;					
			} 

			// if undercompressed, root search up
			if (uc){
				phiL = phi;
				dphi = 0.5*(phiH + phiL) - phi;

				cout << endl;
				cout << "relaxed state found!" << endl;
				cout << "phiL set at nt = " << t << endl;
				cout << "phiH = " << phiH << endl;
				cout << "phiL = " << phiL << endl;
				cout << "dphi = " << dphi << endl;
				cout << "old phi = " << phi << endl;
				this->scale_sys(dphi);
				cout << "new phi = " << phi << endl;
				cout << "increasing root search guess..." << endl;					
			} 

			// if marginal, grow, reset root search
			if (marginal){
				phiL = -1;
				phiH = -1;
				dphi = (1+0.05*drand48())*dphi0;

				cout << endl;
				cout << "marginal state found..." << endl;
				cout << "root search reset at nt = " << t << endl;
				cout << "phiH = " << phi << endl;
				cout << "phiL = " << phiL << endl;
				cout << "dphi = " << dphi << endl;
				cout << "old phi = " << phi << endl;
				this->scale_sys(dphi);
				cout << "new phi = " << phi << endl;
				cout << "marginal, so trying root search again..." << endl;
			}

			// if jammed, end!
			if (jammed){
				cout << endl;
				cout << endl;
				cout << "Found jammed state!" << endl;
				cout << "Writing config at t = " << t*dt << endl;
				cout << "Writing config at nt = " << t << endl;
				cout << "Final U = " << U << endl;
				cout << "Final K = " << K << endl;
				cout << "Final phi = " << phi << endl;
				cout << "Final csum = " << csum << endl;
				cout << "Final niso = " << niso << endl;
				cout << "Final niso max = " << DOF*N-NDIM+1 << endl;
				cout << "Final rattler # = " << nr << endl;
				cout << "Final contacts:" << endl;
				for (int i=1; i<N+1; i++){
					cout << setw(6) << pc[i-1];
					if (i % 10 == 0){
						cout << endl;
					}
				}
				cout << endl;
				if (N < 40){
					cout << "Contact Matrix:" << endl;
					this->print_c_mat();
				}
				if (xyzobj.is_open())
					this->print_xyz();

				cout << endl;
				cout << endl;
				isjammed = 1;
			}
		}
	}

	// test for stalled growth
	if (abs(dphi) < 1e-14){
		phiL = -1;
		phiH = -1;
		dphi = (0.5+drand48())*0.1*dphi0;
		check_rattlers = 0;

		cout << endl;
		cout << "stalled growth found..." << endl;
		cout << "root search reset at nt = " << t << endl;
		cout << "phiH = " << phi << endl;
		cout << "phiL = " << phiL << endl;
		cout << "dphi = " << dphi << endl;
		cout << "old phi = " << phi << endl;
		this->scale_sys(dphi);
		cout << "new phi = " << phi << endl;
		cout << "marginal, so trying root search again..." << endl;
	}
}



void packing::catch_nans(int t){
	// test dynamic variables for nan
	double U0, K0, v0, x0;
	bool test;

	U0 = U;
	K0 = K;
	v0 = v[0][0];
	x0 = x[0][0];
	
	test = std::isnan(U0);
	test = test || std::isnan(K0);
	test = test || std::isnan(v0);
	test = test || std::isnan(x0);

	
	if (test){
		cout << "\n\n\n&&&&& FOUND NANS at t = " << t << " &&&&&&\n\n\n";
		this->print_vars();
		throw "nan found in dynamic variable, ending...\n";
	}
}


/* 

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	FIRE (fast inertial relaxation engine) is implemented instead of damping
	Algorithm is very simple and easy to implement. For more info, see the following paper:

	Bitzek, Erik, Pekka Koskinen, Franz Gähler, Michael Moseler, and Peter Gumbsch. 
	“Structural Relaxation Made Simple.” Physical Review Letters 97, no. 17 (October 27, 2006): 170201.
	DOI: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

*/


void packing::fire(){

	int i,j,d;
	double P = 0;
	double vstarnrm = 0;
	double fstarnrm = 0;
	int Nmin = 50;


	// calculate P
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			P += F[i][d]*v[i][d];
			vstarnrm += v[i][d]*v[i][d];
			fstarnrm += F[i][d]*F[i][d];
		}		
	}

	vstarnrm = sqrt(vstarnrm);
	fstarnrm = sqrt(fstarnrm);

	// update v if forces acting
	if (fstarnrm > 0){
		for (i=0; i<N; i++){
			for (d=0; d<NDIM; d++)
				v[i][d] = (1-alpha)*v[i][d] + alpha*(F[i][d]/fstarnrm)*vstarnrm;
		}
	}

	// now update alphas for P
	if (P >= 0 && np > Nmin){
		// increase dt
		if (dt*finc < dtmax)
			dt *= finc;
		else
			dt = dtmax;

		// decrease alpha
		alpha *= falpha;

		np++;
	}		
	else if (P < 0){
		// reset K to measure based on new info
		K = 0;

		// decrease time step
		dt *= fdec;
		
		// set global velocity vector to zero
		for (i=0; i<N; i++){
			for (d=0; d<NDIM; d++){
				K += 0.5*m[i]*v[i][d]*v[i][d];
				v[i][d] = 0;
			}
		}		

		// set alpha -> alphaStart
		alpha = alpha0;

		// set np -> 0
		np = 0;
	}
	else if (P >= 0 && np <= Nmin)
		np++;


}

/* 
==================================

		GROWTH FUNCTIONS		 

================================== 
*/

void packing::scale_sys(double dphi){
	int i,d;
	double g,msum,vol;	

	// get scale factor
	g = pow((phi+dphi)/phi,(double)1/NDIM);

	// increase phi
	// phi = phi + dphi;

	// scale energy
	ep *= g*g;	

	// scale lengths
	msum = 0;
	for (i=0; i<N; i++){
		// scale particle sizes
		m[i] *= pow(g,NDIM);
		r[i] *= g;
		msum += m[i];
	}

	vol = 1;
	for (d=0; d<NDIM; d++)
		vol *= L[d];
	phi = msum/vol;

	dt *= pow(g,NDIM/2);

	// scale cutoff length
	rcut *= g;	
}





/* 
=========================================

	CELL LIST/NEIGHBOR LIST FUNCTIONS		 

========================================= 
*/

void packing::initialize_nlcl(){
	cout << "** INIALIZING NLCL:" << endl;
	cout << "** getting cell neighbors..." << endl;
	this->get_cell_neighbors();

	cout << "** getting cell positions..." << endl;
	this->get_cell_positions();
	this->print_cell_pos();

	cout << "** populating cell information..." << endl;
	this->update_cell();
	this->print_cell();
	this->print_clabel();
	this->print_celln();	
	this->print_cell_neighbors();	

	cout << "** populating neighbor list information..." << endl;
	this->update_neighborlist();
	this->print_neighborlist();
	cout << "** NLCL INITIALIZATION COMPLETE." << endl;
}


void packing::add_cell(int m, int i){
	cell[m].push_back(i);
}

void packing::reset_cells(){
	int i;

	for (i=0; i<NCELLS; i++)
		cell[i].clear();
}

void packing::add_neighbor(int i, int j){
	neighborlist[i].push_back(j);
}

void packing::reset_neighborlist(int i){
	neighborlist[i].clear();
}

void packing::update_celln(){
	int m;
	for (m=0; m<NCELLS; m++)
		celln[m] = cell[m].size();
}

void packing::get_cell_neighbors(){
	int i,zi;
	if (NDIM == 3){
		for (i=0; i<NCELLS; i++){
			zi = i/(NCL*NCL);

			// faces
			cellneighbors[i][0] = (i+NCELLS-1) % NCELLS; 			// left neighbor (i-1)
			cellneighbors[i][1] = i-NCL;						// bottom neighbor (j-1)
			cellneighbors[i][2] = (i+1) % NCELLS; 				// right neighbor (i+1)
			cellneighbors[i][3] = i+NCL;						// top neighbor (j+1)
			cellneighbors[i][4] = (i+NCELLS-NCL*NCL) % NCELLS;			// backward neighbor (k-1)
			cellneighbors[i][5] = (i+NCL*NCL) % NCELLS;				// forward neighbor (k+1)		

			// y-direction bc
			if (i-zi*NCL*NCL < NCL)
				cellneighbors[i][1] = i-NCL+NCL*NCL;	
			if ((i+NCL)/(NCL*NCL) > zi)
				cellneighbors[i][3] = i-NCL*NCL+NCL;	

			// edges

			// * xy plane
			cellneighbors[i][6] = cellneighbors[i][1] - 1;	// bottom-left neightbor (i-1,j-1)
			cellneighbors[i][7] = cellneighbors[i][1] + 1;	// bottom-right neighbor (i+1,j-1)
			cellneighbors[i][8] = cellneighbors[i][3] - 1; 	// top-left neighbor (i-1,j+1)
			cellneighbors[i][9] = cellneighbors[i][3] + 1;	// top-right neighbor (i+1.j+1)

			// * xz plane
			cellneighbors[i][10] = (cellneighbors[i][4] - 1) % NCELLS;	// back-left neighbor (i-1,k-1)
			cellneighbors[i][11] = (cellneighbors[i][4] + 1) % NCELLS;	// back-right neighbor (i+1,k-1)
			cellneighbors[i][12] = (cellneighbors[i][5] - 1) % NCELLS;	// front-left neighbor (i-1,k+1)		
			cellneighbors[i][13] = (cellneighbors[i][5] + 1) % NCELLS;	// back-front neighbor (i+1,k+1)

			// * yz plane
			cellneighbors[i][14] = cellneighbors[i][4] - NCL;	// back-bottom neighbor (j-1,k-1)
			cellneighbors[i][15] = cellneighbors[i][4] + NCL;	// back-top neighbor (j+1,k-1)
			cellneighbors[i][16] = cellneighbors[i][5] - NCL;	// front-bottom neighbor (j-1,k+1)
			cellneighbors[i][17] = cellneighbors[i][5] + NCL; 	// front-top neighbor (j+1,k+1)

			// y-direction bc
			if (i-zi*NCL*NCL < NCL){
				cellneighbors[i][14] = cellneighbors[i][4]-NCL+NCL*NCL;
				cellneighbors[i][16] = cellneighbors[i][5]-NCL+NCL*NCL;	
			}
			if ((i+NCL)/(NCL*NCL) > zi){
				cellneighbors[i][15] = cellneighbors[i][4]-NCL*NCL+NCL;
				cellneighbors[i][17] = cellneighbors[i][5]-NCL*NCL+NCL;			
			}
			


			// cubic vertices
			cellneighbors[i][18] = (cellneighbors[i][14] - 1) % NCELLS; // back-bottom-left neighbor (i-1,j-1,k-1)
			cellneighbors[i][19] = (cellneighbors[i][14] + 1) % NCELLS; // back-bottom-right neighbor (i+1,j-1,k-1)
			cellneighbors[i][20] = (cellneighbors[i][15] - 1) % NCELLS; // back-top-left neighbor (i-1,j+1,k-1)
			cellneighbors[i][21] = (cellneighbors[i][15] + 1) % NCELLS; // back-top-right neighbor (i+1,j+1,k-1)
			cellneighbors[i][22] = (cellneighbors[i][16] - 1) % NCELLS; // front-bottom-left neighbor (i-1,j-1,k+1)
			cellneighbors[i][23] = (cellneighbors[i][16] + 1) % NCELLS; // front-bottom-right neighbor (i+1,j-1,k+1)
			cellneighbors[i][24] = (cellneighbors[i][17] - 1) % NCELLS; // front-top-left neighbor (i-1,j+1,k+1)
			cellneighbors[i][25] = (cellneighbors[i][17] + 1) % NCELLS; // front-top-right neighbor (i+1,j+1,k+1)



			if (i % NCL == 0){
				// left BC
				cellneighbors[i][0] = i+NCL-1;
				cellneighbors[i][6] = cellneighbors[i][1]+NCL-1;
				cellneighbors[i][8] = cellneighbors[i][3]+NCL-1;
				cellneighbors[i][10] = cellneighbors[i][4]+NCL-1;
				cellneighbors[i][12] = cellneighbors[i][5]+NCL-1;	
				cellneighbors[i][18] = cellneighbors[i][14]+NCL-1;
				cellneighbors[i][20] = cellneighbors[i][15]+NCL-1;
				cellneighbors[i][22] = cellneighbors[i][16]+NCL-1;
				cellneighbors[i][24] = cellneighbors[i][17]+NCL-1;
			}
			if ((i+1) % NCL == 0){
				// right BC
				cellneighbors[i][2] = i-NCL+1;
				cellneighbors[i][7] = cellneighbors[i][1]-NCL+1;
				cellneighbors[i][9] = cellneighbors[i][3]-NCL+1;
				cellneighbors[i][11] = cellneighbors[i][4]-NCL+1;
				cellneighbors[i][13] = cellneighbors[i][5]-NCL+1;
				cellneighbors[i][19] = cellneighbors[i][14]-NCL+1;
				cellneighbors[i][21] = cellneighbors[i][15]-NCL+1;
				cellneighbors[i][23] = cellneighbors[i][16]-NCL+1;
				cellneighbors[i][25] = cellneighbors[i][17]-NCL+1;
			}
		}
	}	
	else if(NDIM == 2){
		cout << "ERROR: NDIM=2 not supported when using neighborlist..." << endl;
		throw "NDIM not supported!\n";
	}
	else{
		cout << "ERROR: NDIM not supported when using neighborlist..." << endl;
		throw "NDIM not supported!\n";
	}
}

void packing::get_cell_positions(){
	int m,d;
	int ind;	

	// loop over cells, get pos of each one
	for (m=0; m<NCELLS; m++){
		for (d=0; d<NDIM; d++){
			// get indices xi,yi,zi,wi,etc
			ind = floor((m % (int)pow(NCL,d+1))/(pow(NCL,d)));

			// get position in given dimension of cell
			cellpos[m][d] = g[d]*(ind+0.5);
		}
	}
}

int packing::get_new_cell(int i){
	int m,d,insq,minm;
	double cposi;
	double dr,h,mindist;

	mindist = 3*L[0];
	minm = -1;

	// loop over cells, get pos of each one
	for (m=0; m<NCELLS; m++){

		h = 0;
		insq = 0;
		for (d=0; d<NDIM; d++){
			// get cell pos
			cposi = cellpos[m][d];

			// get distance to particle i
			dr = x[i][d] - cposi;
			dr = dr - L[d]*round(dr/L[d]);

			h += dr*dr;					
		}
		h = sqrt(h);

		// check if min dist
		if (h < mindist){			
			mindist = h;
			minm = m;
		}
	}	

	// if no min set, error
	if (minm == -1){
		cout << "ERROR: No native cell found, error in get_new_cell()...";
		throw "NO NATIVE CELL FOUND!!\n";
	}
	else
		return minm;
}

void packing::update_cell(){
	// local variables
	int i,d,newcell;

	// reset cells
	this->reset_cells();

	// for all atoms in system, check cell location
	for (i=0; i<N; i++){
		// get updated location of particle		
		newcell = this->get_new_cell(i);		

		// append particle to new cell
		this->add_cell(newcell,i);
		clabel[i] = newcell;		
	}

	// get new cell numbers
	this->update_celln();
}


void packing::update_neighborlist(){
	// local variables
	int M,i,j,d,l,l2,m;
	int p2;
	double dr,h;

	// loop over atoms in system
	for (i=0; i<N; i++){		

		// clear neighborlist, to be updated
		this->reset_neighborlist(i);

		// get cell loc of atom i
		l = clabel[i];

		// get number of atoms in cell l
		M = celln[l];

		// cout << "checking NL for i = " << i;
		// cout << "; in cell " << l;
		// cout << "; # of atoms = " << M;	
		// cout << endl;	

		// loop over atoms in cell l
		for (j=0; j<M; j++){

			// get particle in cell 
			p2 = cell[l].at(j);

			// if diff from particle i, check distance
			if (p2 != i){
				h = this->get_distance(i,p2);
				if (h < rcut){
					// cout << "found neighbor of i = " << i << ", at p2 = " << p2 << endl;
					this->add_neighbor(i,p2);
				}
			}
		}

		// cout << "; checking neighboring cells";

		// loop over cells neighboring to cell l
		for (m=0; m<NCN; m++){
			// get cell label of neighboring cell to cell l
			l2 = cellneighbors[l][m];

			// get number of atoms in cell l2
			M = celln[l2];

			// cout << "; l2 = " << l2;
			// cout << ", M = " << M;

			// loop over atoms in cell l2
			for (j=0; j<M; j++){

				// get particle in cell
				p2 = cell[l2].at(j);

				// if diff from particle i, check distance
				if (p2 != i){
					h = this->get_distance(i,p2);
					if (h < rcut)
						this->add_neighbor(i,p2);
				}
			}
		}
	}
}































