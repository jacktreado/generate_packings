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


/*	FORCE UPDATE  */

double packing::hs(double sij, double xij){
	double f;

	// calc pairwise force, update potential
	f = (ep/(sij*sij))*(1-sij/xij);
	// f = ep*(1-sij/xij);

	return f;
}

void packing::hs_force(){
	int i,j,d,cind,jstart,jend,jj;
	double sij, dx;
	double xij[NDIM];
	bool use_ncnl;

	// check whether or not to use neighbor list
	use_ncnl = (NCL > -1);
	if (use_ncnl)
		jstart = 0;
	else
		jend = N;

	// reset U and contacts
	U = 0;
	this->reset_c();


	for (i=0; i<N; i++){
		if (use_ncnl)
			jend = neighborlist[i].size();
		else
			jstart = i+1;

		for (jj=jstart; jj<jend; jj++){
			// get particle index if using nlcl
			if (use_ncnl)
				j = neighborlist[i].at(jj);
			else
				j = jj;

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
		}
	}
}


void packing::vel_update(){
	int i,d;
	double vel, anew;

	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			vel = v[i][d];

			// update velocity
			anew = F[i][d]/m[i];
			vel += 0.5*dt*(anew + aold[i][d]);
			v[i][d] = vel;
			aold[i][d] = anew;
		}
	}
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
		// count contacts
		r = pc[i];

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
	tmp0 = 0.01;

	// initialize particle velocities
	cout << "* initializing velocities in jamming_finder()..." << endl;
	this->rand_vel_init(tmp0);

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
	dUtol = 1e-8;
	epc = 0;
	epcN = 5e2;

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
		// update nearest neighbor lists if applicable
		if (t % nnupdate == 0 && NCL > -1){
			cout << "^ ";
			this->update_nlcl(t);
		}

		// first md time step
		this->pos_update();

		// force update
		this->hs_force();

		// fire update (inertial relaxer)
		this->fire();

		// vel update
		this->vel_update();		

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
			if (xyzobj.is_open()){				
				if (NCL > -1)
					this->print_nl_xyz();
				else
					this->print_xyz();				
			}
			if (enobj.is_open())
				this->print_en(t);
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
	bool gr,oc,uc,jammed;
	double dphi = dphi0;

	gr = 0;
	oc = 0;
	uc = 0;
	jammed = 0;

	nbb = N - nr;
	niso = DOF*nbb - NDIM + 1;
	csum = this->get_c_sum();

	gr = (U < Utol);
	oc = (U > 2 * Utol && K < Ktol && epconst == 1 && csum > 0);
	uc = (U < Utol && epconst == 1);
	jammed = ( (U > Utol && U < 2 * Utol && K < Ktol && epconst == 1 && csum > 0) );


	if (phiH < 0){
		if (gr){
			check_rattlers = 0;
			this->scale_sys(dphi);
		}
		else if (oc){				
			phiH = phi;
			dphi = -dphi0;
			check_rattlers = 1;

			cout << endl;
			cout << "phiH 1st set at nt = " << t << endl;
			this->monitor_scale(dphi,phiH,phiL);
			cout << "entering root search..." << endl;
			cout << endl;

			// if NLCL, change update check (particles don't move, don't need to check as often)
			if (NCL > -1)
				nnupdate *= 50;
		}
	}
	else{
		if (phiL < 0){

			// if still overcompressed, decrease again
			if (oc){				
				phiH = phi;
				dphi = -drand48() * dphi0;

				cout << endl;
				cout << "still overcompressed..." << endl;
				cout << "phiH set at nt = " << t << endl;
				this->monitor_scale(dphi,phiH,phiL);				
				cout << "continuing root search..." << endl;
			}

			// if undercompressed, set phiL, root search
			if (uc){
				phiL = phi;
				dphi = 0.5*(phiH + phiL) - phi;

				cout << endl;
				cout << "relaxation found!" << endl;
				cout << "phiL set at nt = " << t << endl;
				this->monitor_scale(dphi,phiH,phiL);
				cout << "root search interval found, growing..." << endl;
			}

			if (jammed){
				phiL = 0.99*phi;
				dphi = 0.5*(phiH + phiL) - phi;

				cout << endl;
				cout << "almost jammed found!" << endl;
				cout << "phiL set at nt = " << t << endl;
				this->monitor_scale(dphi,phiH,phiL);
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
				this->monitor_scale(dphi,phiH,phiL);
				cout << "decreasing root search guess..." << endl;					
			} 

			// if undercompressed, root search up
			if (uc){
				phiL = phi;
				dphi = 0.5*(phiH + phiL) - phi;

				cout << endl;
				cout << "relaxed state found!" << endl;
				cout << "phiL set at nt = " << t << endl;
				this->monitor_scale(dphi,phiH,phiL);
				cout << "increasing root search guess..." << endl;					
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
				cout << "Final phi = " << setprecision(10) << phi << endl;
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

				if (xyzobj.is_open()){
					if (NCL < 0)		
						this->print_xyz();
					else{
						this->print_nl_xyz();
						this->print_all_nl_xyz();
					}

				}
				cout << endl;
				cout << endl;
				isjammed = 1;
			}
		}
	}

	// // test for stalled growth
	// if (abs(dphi) < 1e-14){
	// 	phiL = -1;
	// 	phiH = -1;
	// 	dphi = 0.05 * (2 * drand48() - 1) * dphi0;
	// 	check_rattlers = 0;

	// 	cout << endl;
	// 	cout << "stalled growth found..." << endl;
	// 	cout << "root search reset at nt = " << t << endl;
	// 	this->monitor_scale(dphi,phiH,phiL);
	// 	cout << "marginal, so trying root search again..." << endl;
	// }
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
		cout << "* ";
		// reset K to measure based on new info
		K = 0;

		// decrease time step
		if (dt * fdec > 1e-2 * dtmax)
			dt *= fdec;
		else
			dt = 1e-2 * dtmax;
		
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
	double s,msum,vol;	

	// get scale factor
	s = pow((phi+dphi)/phi,(double)1/NDIM);

	// scale lengths
	msum = 0;
	for (i=0; i<N; i++){
		// scale particle sizes
		m[i] *= pow(s,NDIM);
		r[i] *= s;
		msum += m[i];
	}

	vol = 1;
	for (d=0; d<NDIM; d++)
		vol *= L[d];
	phi = msum/vol;

	ep *= pow(s, 2);
	dt *= pow(s, 0.5 * NDIM);
	dtmax *= pow(s, 0.5 * NDIM);

	// if NLCL engaged, scale rcut
	if (NCL > -1) {
		this->scale_rcut(s);
		if (NCL > 3)
			this->update_cell();

		this->update_neighborlist();
	}
}





/* 
=========================================

	CELL LIST/NEIGHBOR LIST FUNCTIONS		 

========================================= 
*/

void packing::init_rcut(){
	// check if rcut pts to null
	if (rcut == nullptr){
		cout << "rcut not yet initialized, initializing to N..." << endl;
		rcut = new double[N];
	}

	int i;
	double maxrad;
	maxrad = 0;

	for (i=0; i<N; i++){
		if (r[i] > maxrad)
			maxrad = r[i];
	}

	// loop over system, make rcut avg between maxrad and ri
	for (i=0; i<N; i++)
		rcut[i] = r[i] + maxrad;

}

void packing::scale_rcut(double s){
	int i;
	for (i=0; i<N; i++)
		rcut[i] *= s;
}

void packing::update_nlcl(int t){
	if (NCL > 3)
		this->update_cell();

	this->update_neighborlist();
	this->nlcl_check(t);
}

void packing::add_cell(int k, int i){
	cell[k].push_back(i);
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
	int k;
	for (k=0; k<NCELLS; k++)
		celln[k] = cell[k].size();
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
	int k,d,i1,i2;
	int ind;

	// loop over cells, get pos of each one
	for (k=0; k<NCELLS; k++){		
		for (d=0; d<NDIM; d++){			
			// get indices xi,yi,zi,wi,etc
			i1 = round(pow(NCL,d+1));
			i2 = round(pow(NCL,d));
			ind = floor((k % i1)/(i2));

			// get position in given dimension of cell
			cellpos[k][d] = g[d]*(ind+0.5);
		}
	}
}

int packing::get_new_cell(int i){
	int k,d,insq,minm;
	double cposi;
	double dr,h,mindist;

	mindist = 3*L[0];
	minm = -1;

	// loop over cells, get pos of each one
	for (k=0; k<NCELLS; k++){

		h = 0;
		insq = 0;
		for (d=0; d<NDIM; d++){
			// get cell pos
			cposi = cellpos[k][d];

			// get distance to particle i
			dr = x[i][d] - cposi;
			dr = dr - L[d]*round(dr/L[d]);

			h += dr*dr;					
		}
		h = sqrt(h);

		// check if min dist
		if (h < mindist){
			mindist = h;
			minm = k;			
		}
	}

	// if no min set, error
	if (minm == -1){
		cout << "ERROR: No native cell found, error in get_new_cell()...";
		cout << "ERROR: particle position = " << x[i][0] << ", " << x[i][1] << ", " << x[i][2] << endl;
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
	int K,i,j,d,l,l2,k;
	int p2;
	double dr,h;

	// loop over atoms in system
	for (i=0; i<N; i++){		

		// clear neighborlist, to be updated
		this->reset_neighborlist(i);

		// get cell loc of atom i
		l = clabel[i];

		// get number of atoms in cell l
		K = celln[l];	

		// loop over atoms in cell l
		for (j=0; j<K; j++){

			// get particle in cell 
			p2 = cell[l].at(j);

			// if p2 > particle i, check distance (knocks down force update loop)
			if (p2 > i){
				h = this->get_distance(i,p2);
				if (h < rcut[i])
					this->add_neighbor(i,p2);
			}
		}

		// loop over cells neighboring to cell l
		for (k=0; k<NCN; k++){
			// get cell label of neighboring cell to cell l
			l2 = cellneighbors[l][k];

			// get number of atoms in cell l2
			K = celln[l2];			

			// loop over atoms in cell l2
			for (j=0; j<K; j++){

				// get particle in cell
				p2 = cell[l2].at(j);

				// if p2 > particle i, check distance (knocks down force update loop)
				if (p2 > i){
					h = this->get_distance(i,p2);
					if (h < rcut[i])
						this->add_neighbor(i,p2);
				}
			}
		}
	}
}



void packing::nlcl_check(int t){
	int i,j,k,nf;
	double sij,h;

	for (i=0; i<N; i++){
		for (j=i+1; j<N; j++){
			// get contact distance sij
			sij = r[j] + r[i];

			// get distance between particles
			h = this->get_distance(i, j);

			// if distance < sij, check that two particles are neighbors
			nf = 0;
			if (h < sij){

				// check neighbor list of i
				K = neighborlist[i].size();
				for (k=0; k<K; k++){
					if (neighborlist[i][k]==j){
						nf = 1;
						break;
					}
				}

				// if j not in nl of i, then problem!
				if (nf == 0){
					this->print_cell_neighbors();
					if (xyzobj.is_open())
						this->print_nl_xyz(i,j);
					cout << "*** t = " << t << ", neighbor not picked up in nlcl:" << endl;
					cout << "i = " << i << ", j = " << j << ", dist/rcut[i] = " << h/rcut[i] << endl;
					cout << "cells: ci = " << clabel[i] << ", cj = " << clabel[j] << endl;
					cout << "ending program..." << endl;
					throw "neighbor not caught\n";
				}
			}
		}
	}
}



























