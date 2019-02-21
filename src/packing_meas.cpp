/*

	Methods implementation 
	for packing class

	BY Jack Treado

*/

#include "packing.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
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

// FUNCTION to populate entries in dynamical matrix
// assuming that forces are already calculated

/*
void packing::populate_dynamical_matrix(double h0, string& dmstr){
	// local variables
	int i,j,d,k,l,e,d1,d2,h;
	double Mkl;
	h = h0/L;

	int Nentries = (NDIM*N*(NDIM*N+1))/2; 				// number of unique entries in dynamical matrix
	vector<double> Fk(NDIM*N,0);						// vector of all unperturbed forces
	Eigen::MatrixXd P(NDIM*N,NDIM*N);			// matrix of all perturbed forces	
	vector<double> M(Nentries,1);						// dynamical matrix M

	// store all unperturbed forces first
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			k = NDIM*i + d;
			Fk.at(k) = F[i][d];
		}
	}

	// get all perturbed forces
	for (j=0; j<N; j++){

		// loop over perturbation directions
		for (d2=0; d2<NDIM; d2++){

			// get linear index
			l = NDIM*j+d2;

			// perturb particle j in d2 direction, get updated forces
			this->perturbed_force(j,d2,h);

			// loop over to extract forces to perturb
			for (i=0; i<N; i++){
				// perturb in NDIM dimensions
				for (d1=0; d1<NDIM; d1++){
					// get linear index
					k = NDIM*i+d1;

					// add forces to matrix
					P(k,l) = F[i][d1];
				}
			}

		}
	}


	// loop over particles, get forces as a function of particle perturbations
	// definition from: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
	for (i=0; i<N; i++){

		// loop over particles j
		for (j=i; i<N; j++){

			// loop over dimensions
			for (d1=1; d1<NDIM; d1++){
				for (d2=1; d2<NDIM; d2++){

					// entry in force matrix Pkl
					k = NDIM*i + d1;
					l = NDIM*j + d2;

					// calculate entry in dynamical matrix
					if (k > l)
						e = (3*NDIM)*l + k - ((l+1)*l)/2;
					else
						e = (3*NDIM)*k + l - ((k+1)*k)/2;

					// calculate single matrix entry
					Mkl = (Fk(k) - P(k,l) + Fk(l) - P(l,k))/(2*h);

					// update matrix entry
					M.at(e) = P(k,l);
				}
			}
		}
	}

	// print dynamical matrix to file

}


void packing::perturbed_force(int i, int d, double h){
	// local variables
	double V,x0;

	// perturb particle
	x0 = x[i][d];
	x[i][d] += h;

	// calc new force
	this->hs_force();

	// reverse perturbation
	x[i][d] = x0;
}


*/



























void packing::calc_dm(string& fstr){
	// local variables
	int i,k,j,l,d,nr,kr;
	int p,q;
	int df = DOF*N;
	double** dm;
	double V,Vi,Vj,Vij,Hij;
	const double h = 1e-8;	

	// intialize dynamical matrix
	dm = new double*[df];
	for (p=0; p<df; p++)
		dm[p] = new double[df];

	// get unperturbed V
	U = 0;
	this->hs_force();

	// get number of rattlers
	kr = 0;
	nr = this->rmv_rattlers(kr);

	// initialize perturbed V's
	V = U;
	Vi = 0;
	Vj = 0;
	Vij = 0;
	
	// loop over particles	
	for (i=0; i<N; i++){

		// perturb particle i
		for (k=0; k<DOF; k++){
			Vi = this->perturb_single_particle(i,k,h);

			// get matrix index
			p = i*DOF+k;

			// loop over particles j, perturb dof l
			for (j=i; j<N; j++){

				// perturb particle j
				for (l=0; l<DOF; l++){
					// get matrix index
					q = j*DOF+l;

					if (p > DOF*N || q > DOF*N){
						cout << "over accessing dm matrix!" << endl;
						throw "out of bounds\n";
					}

					// get perturbed particle
					Vj = this->perturb_single_particle(j,l,h);
					Vij = this->perturb_two_particles(i,j,k,l,h);
					Hij = (Vij-Vi-Vj+V)/(h*h);

					// save to dynamical matrix
					dm[p][q] = Hij;
					dm[q][p] = Hij;
				}
			}
		}
	}

	// print dm to file, free up memory
	int w = 32;
	int pr = 10;	
	ofstream obj(fstr.c_str());
	obj << df << endl;
	obj << nr << endl;
	for(p=0; p<df; p++){
		for (q=0; q<df; q++)
			obj << setprecision(pr) << setw(w) << dm[p][q];
		obj << endl;
		delete [] dm[p];
	}
	obj.close();
	delete [] dm;
}

double packing::perturb_single_particle(int i, int k, double h){
	double V,x0;

	// perturb particle
	x0 = x[i][k];
	x[i][k] += h;

	// calc new system
	U = 0;
	this->hs_force();

	// get energy
	V = U;

	// reverse perturbation
	x[i][k] = x0;

	// return V
	return V;
}

double packing::perturb_two_particles(int i, int j, int k, int l, double h){
	double V,x10,x20;	

	// perturb particle
	x10 = x[i][k];
	x20 = x[j][l];
	x[i][k] += h;
	x[j][l] += h;

	// calc new system
	U = 0;
	this->hs_force();

	// get energy
	V = U;

	// reverse perturbation
	x[i][k] = x10;
	x[j][l] = x20;

	// return V
	return V;
}

// calculate vacf from already initialized system
void packing::calc_vacf(int NT, int vsave, ofstream& obj){
	// get number of time steps
	int i,d,t;

	// initialize array to save velocities
	vector<double>** vlist;
	vlist = new vector<double>*[N];
	for (i=0; i<N; i++)
		vlist[i] = new vector<double>[NDIM];

	// initialize vacf array
	int vsamp = NT/vsave;
	vector<double> vacf;
	vector<double> numer;
	vector<double> denom;
	for (i=0; i<vsamp; i++){
		numer.push_back(0.0);
		denom.push_back(0.0);
		vacf.push_back(0.0);
	}

	// loop over time
	cout << "* looping over time, saving velocities" << endl;
	for (t=0; t<NT; t++){
		// first md step
		this->pos_update();

		// force update
		this->hs_force();

		// vel update
		this->vel_update();

		// get vacf at right time
		if (t % vsave == 0){
			this->get_vacf(t,vsave,vlist,numer,denom);			
		}

		// print info
		if (t % plotskip == 0){
			this->md_monitor(t,-1,-1,-1);
			cout << endl << endl;
			if (xyzobj.is_open())
				this->print_xyz();
		}
	}

	this->finish_vacf(NT,vsave,numer,denom,vacf);

	// delete vlist (velocity saver)
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			vlist[i][d].clear();

		delete [] vlist[i];
	}
	delete [] vlist;

	cout << "VACF MD completed!" << endl;
	this->print_vacf(NT,vsave,vacf,obj);
}

// get vacf from previous velocities
void packing::get_vacf(int t, int vsave, vector<double>** vlist, vector<double>& numer, vector<double>& denom){	
	int i,j,d,dt,dtmaxtmp;

	// push back current velocities
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			vlist[i][d].push_back(v[i][d]);
	}

	// loop over previous velocities, calc vacf
	dt = 1;
	dtmax = t/vsave+1;	
	while (dt < dtmax){		
		for (i=0; i<N; i++){
			for (d=0; d<NDIM; d++){
				numer.at(dt-1) += (vlist[i][d].at(dtmax-1)*vlist[i][d].at(dtmax-dt-1));
				denom.at(dt-1) += (vlist[i][d].at(dtmax-dt-1)*vlist[i][d].at(dtmax-dt-1));
			}
		}

		// increase dt	
		dt++;
	}
}

void packing::finish_vacf(int NT, int vsave, vector<double>& numer, vector<double>& denom, vector<double>& vacf){
	int i,j,d,dt,dtmaxtmp;

	// loop over previous velocities, calc vacf
	dt = 1;
	dtmax = NT/vsave;	

	while (dt < dtmax){				
		// get quotient
		vacf.at(dt-1) = numer.at(dt-1)/denom.at(dt-1);

		// increase dt	
		dt++;			
	}
}

void packing::print_vacf(int NT, int vsave, vector<double>& vacf, ofstream& obj){
	int s = vacf.size();
	int i,nsamp;
	nsamp = NT/vsave;

	obj << nsamp << endl;
	obj << dt*vsave << endl;
	obj << dt << endl;
	for (i=0; i<s; i++)
		obj << vacf.at(i) << endl;
}








