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


void packing::dynamical_matrix(string& dmstr, double h0){
	// local variables
	int i,j,d,k,l,e,d1,d2,nr,kr;
	double Mkl,h,dt0,len;

	// rescale all lengths to smaller radius
	len = r[0];
	this->rescale_lengths(len);

	// set h scale
	h = h0*L[0];

	// remove any rattlers
	kr = 0;
	nr = this->rmv_rattlers(kr);
	cout << "there are nr = " << nr << " rattlers" << endl;

	int Nentries = (NDIM*N*((NDIM*N)+1))/2; 			// number of unique entries in dynamical matrix
	vector<double> Fk(NDIM*N,0);						// vector of all unperturbed forces
	Eigen::MatrixXd Pplus(NDIM*N,NDIM*N);				// matrix of all perturbed force (+h)
	Eigen::MatrixXd Pminus(NDIM*N,NDIM*N);				// matrix of all perturbed force (-h)
	Eigen::MatrixXd numHessian(NDIM*N,NDIM*N);			// numerical hessian
	Eigen::MatrixXd massMatrix(NDIM*N,NDIM*N);			// mass matrix
	Eigen::MatrixXd anHessian(NDIM*N,NDIM*N);			// analytical hessian


	// store all unperturbed forces first
	cout << "getting forces..." << endl;	
	this->hs_force();
	cout << "total potential energy U = " << U << endl;
	for (i=0; i<N; i++){		

		cout << "Force on particle i = " << i << ": ";
		for (d=0; d<NDIM; d++){
			k = NDIM*i + d;

			// only count force if not a rattler	
			// if (pc[i] > 0)		
			// 	Fk.at(k) = F[i][d];
			// else
			// 	Fk.at(k) = 0;
			Fk.at(k) = F[i][d];

			cout << setw(20) << setprecision(6) << Fk.at(k);
		}
		cout << endl;
	}

	// get all perturbed forces
	for (j=0; j<N; j++){

		// loop over perturbation directions
		for (d2=0; d2<NDIM; d2++){

			// get linear index
			l = NDIM*j+d2;

			// perturb particle j in d2 (+h) direction, get updated forces
			this->perturbed_force(j,d2,h);

			// loop over to extract forces to perturb
			for (i=0; i<N; i++){

				// perturb in NDIM dimensions
				for (d1=0; d1<NDIM; d1++){
					// get linear index
					k = NDIM*i+d1;

					// add forces to matrix
					Pplus(k,l) = F[i][d1];
				}
			}

			// perturb particle j in d2 (+h) direction, get updated forces
			this->perturbed_force(j,d2,-h);

			// loop over to extract forces to perturb
			for (i=0; i<N; i++){

				// perturb in NDIM dimensions
				for (d1=0; d1<NDIM; d1++){
					// get linear index
					k = NDIM*i+d1;

					// add forces to matrix
					Pminus(k,l) = F[i][d1];
				}
			}

		}
	}

	// store all unperturbed forces first
	cout << "outputting forces again to check if perturbations remain..." << endl;
	U = 0;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			F[i][d] = 0;
		}
	}
	cout << "reprinting forces after perturbation..." << endl;	
	this->hs_force();
	cout << "total potential energy U = " << U << endl;
	for (i=0; i<N; i++){		

		cout << "Force on particle i = " << i << ": ";
		for (d=0; d<NDIM; d++){
			k = NDIM*i + d;

			// only count force if not a rattler	
			// if (pc[i] > 0)		
			// 	Fk.at(k) = F[i][d];
			// else
			// 	Fk.at(k) = 0;
			// Fk.at(k) = F[i][d];

			cout << setw(20) << setprecision(6) << F[i][d];
		}
		cout << endl;
	}


	// loop over particles, get forces as a function of particle perturbations
	// definition from: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm

	// print dynamical matrix to file
	ofstream obj(dmstr.c_str());
	if (!obj.is_open()){
		cout << "could not open file " << dmstr << endl;
		throw;
	}

	// print out matrix
	obj << N << endl;
	obj << NDIM << endl;
	obj << nr << endl;
	for (k=0; k<NDIM*N; k++){
		for (l=k; l<NDIM*N; l++){
			// get matrix entry
			e = (N*NDIM)*k + l - (((k+1)*k)/2);

			// calculate single matrix entry
			// Mkl = (Fk.at(k) - P(k,l) + Fk.at(l) - P(l,k))/(2*h);
			Mkl = (Pminus(k,l) - Pplus(k,l) + Pminus(l,k) - Pplus(l,k))/((4.0/5.0)*h);

			// store matrix element in dynamical matrix
			numHessian(k,l) = Mkl;
			numHessian(l,k) = Mkl;
		}

		// also populate mass matrix
		i = k % NDIM;
		massMatrix(k,k) = m[i];
	}

	// solve generalized eigenvalue problem for numerical hessian
	cout << "doing eigenvalue decomposition for numerical hessian" << endl;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> num_eigs(numHessian, massMatrix);
    cout << "Generalized eigenvalues for numerical hessian (verify):\n" << num_eigs.eigenvalues() << endl;

    cout << "compare to analytical hessian:" << endl;
    compute_analytical_hessian(anHessian);
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> an_eigs(anHessian, massMatrix);
	cout << "Generalized eigenvalues for analytical hessian:\n" << an_eigs.eigenvalues() << endl;

	// print everything to object
	obj << num_eigs.eigenvalues() << endl;
	obj << an_eigs.eigenvalues() << endl;

	// close object
	obj.close();
}


void packing::perturbed_force(int i, int d, double h){
	// local variables
	int i0,d0;
	double x0;

	// perturb particle
	x0 = x[i][d];
	x[i][d] += h;

	// reset force
	for (i0=0; i0<N; i0++){
		for (d0=0; d0<NDIM; d0++)
			F[i0][d0] = 0.0;
	}

	// calc new force
	this->hs_force();

	// reverse perturbation
	x[i][d] = x0;
}



// analytical method to calculate the DM, to check mixed numerical method
void packing::compute_analytical_hessian(Eigen::MatrixXd& hessian){
	// local variables
	int i,j,jj,l,k,e,d1,d2,Nentries,DF;
	double sij,hij,K1,K2;
	double rij[NDIM];
	double Mkl = 0;

	// get number of total entries
	DF = N*NDIM;
	Nentries = (DF*(DF+1))/2;

	// loop over pairwise particles & directions
	for (k=0; k<DF; k++){
		// get particle i
		i = floor(k/NDIM);

		// get dimension d1
		d1 = k % NDIM;
		for (l=k; l<DF; l++){
			// get particle j
			j = floor(l/NDIM);

			// get dimension d2
			d2 = l % NDIM;

			// get matrix entry
			e = DF*k + l - (((k+1)*k)/2);

			// get matrix element
			Mkl = 0.0;

			// IF i == j, then same particle as perturbation, so calculate every 
			// force due to every contacting particle
			if (i == j){
				for (jj=0; jj<N; jj++){
					if (jj == i)
						continue;

					// get contact distance
					sij = r[i] + r[jj];

					// get real distance
					hij = this->get_distance(i,jj,rij);

					if (hij < sij){
						K1 = ep/(hij*sij);
						K2 = rij[d2]*rij[d1];
						if (d2 == d1)
							K2 += (1 - (hij/sij));
						Mkl += K1*K2;
					}
				}
			}

			// else, only calculate every pairwise force
			else{
				// get contact distance
				sij = r[i] + r[j];

				// get real distance
				hij = this->get_distance(i,j,rij);

				if (hij < sij){
					K1 = ep/(hij*sij);
					K2 = rij[d2]*rij[d1];
					if (d2 == d1)
						K2 += (1 - (hij/sij));
					Mkl = -1*K1*K2;
				}
			}

			// store matrix element in analytical hessian
			hessian(k,l) = Mkl;
			hessian(l,k) = Mkl;
		}
	}
}

// calculate vacf from already initialized system
void packing::calc_vacf(int NT, int vsave, double T0, ofstream& obj){
	// get number of time steps
	int i,d,t;

	// give initial velocities
	this->rand_vel_init(T0);
	
	// initialize vacf arrays
	int vsamp = NT/vsave;
	vector<double> vacf;
	vector<double> numer;
	double denom = 0.0;
	for (i=0; i<vsamp; i++){
		numer.push_back(0.0);
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

		// print info
		if (t % plotskip == 0){
			this->md_monitor(t,-1,-1,-1);
			cout << endl << endl;
			if (xyzobj.is_open())
				this->print_xyz();
		}

		// get vacf at right time
		if (t % vsave == 0){
			this->get_vacf(t,vsave,numer,denom);
		}
	}

	// finish off the vacf!
	this->finish_vacf(NT,vsave,numer,denom,vacf);

	cout << "VACF MD completed!" << endl;
	this->print_vacf(NT,vsave,vacf,obj);

	// close object
	obj.close();
}

// get vacf from previous velocities
void packing::get_vacf(int t, int vsave, vector<double>& numer, double& denom){	
	int i,j,d,tcurrent,t0,dt;
	double v0,vplusdt;

	// push back current velocities
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			vlist[i][d].push_back(v[i][d]);
		}
	}

	// loop over previous velocities, calc vacf
	dt = 0;
	tcurrent = t/vsave;
	while (dt < tcurrent){
		// get dt, time between two points
		t0 = tcurrent - dt;

		// add to numerator
		for (i=0; i<N; i++){
			for (d=0; d<NDIM; d++){
				// velocity at current time t
				v0 = vlist[i][d].at(t0);

				// velocity at time t + dt
				vplusdt = vlist[i][d].at(tcurrent);

				// add to numerator
				numer.at(dt) += v0*vplusdt;
			}
		}

		// increment time 
		dt++;
	}

	// add to denominator
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			denom += v[i][d]*v[i][d];
	}	
}

void packing::finish_vacf(int NT, int vsave, vector<double>& numer, double& denom, vector<double>& vacf){
	int i,j,d,dt,vsamp;

	// average denominator
	denom /= NT*N;

	// loop over previous velocities, calc vacf
	dt = 0;
	vsamp = NT/vsave;	

	while (dt < vsamp){
		// average denominator
		numer.at(dt) /= NT*N;

		// get quotient
		vacf.at(dt) = numer.at(dt)/denom;

		// increase dt	
		dt++;			
	}
}

void packing::print_vacf(int NT, int vsave, vector<double>& vacf, ofstream& obj){
	int s = vacf.size();
	int i,nsamp;
	nsamp = NT/vsave;

	obj << nsamp << endl;
	obj << dt*nsamp << endl;
	obj << dt << endl;
	for (i=0; i<s; i++)
		obj << setw(30) << setprecision(16) << vacf.at(i) << endl;
}




















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








