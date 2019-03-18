/*

	Measure methods implementation
	for rigidbody class

	BY Jack Treado
	10/01/2018

*/

#include "rigidbody.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <vector>

using namespace Eigen;
using namespace std;

// PI
const double PI = 4*atan(1);


void rigidbody::rb_dynamical_matrix(string& dmstr, double h0){
	// local variables
	int i,j,d,k,l,e,d1,d2,nr,kr;
	double Mkl,h,Sij,Iij,am,dnorm,rmin;	

	// rescale system to H atom = 1.0 A, unit energy scale
	rmin = 1e6;
	for (i=0; i<Na[0]; i++){
		if (ar[0][i] < rmin)
			rmin = ar[0][i];
	}
	this->rb_rescale(rmin);
	this->set_ep(10.0);
	h = h0*L[0];

	// remove any rattlers
	kr = 0; 
	nr = this->rmv_rattlers(kr);
	cout << "there are nr = " << nr << " rattlers" << endl;

	int Nentries = (DOF*N*((DOF*N)+1))/2; 			// number of unique entries in dynamical matrix
	Eigen::MatrixXd Pplus(DOF*N,DOF*N);				// matrix of all perturbed force (+h)
	Eigen::MatrixXd Pminus(DOF*N,DOF*N);			// matrix of all perturbed force (-h)

	// matrices to do eigenvalue decomposition in c++
	Eigen::MatrixXd dynMatrix(DOF*N,DOF*N);
	Eigen::MatrixXd massMatrix(DOF*N,DOF*N);	

	// output all forces to check
	cout << "outputting forces..." << endl;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			F[i][d] = 0.0;
			TqW[i][d] = 0.0;
		}
	}
	this->force_update();
	cout << "total potential energy U = " << U << endl;
	for (i=0; i<N; i++){		
		cout << "Force/Torque on particle i = " << i << ": ";
		for (d=0; d<NDIM; d++){
			cout << setw(20) << setprecision(6) << F[i][d];
		}
		for (d=0; d<NDIM; d++){
			cout << setw(20) << setprecision(6) << TqW[i][d];
		}
		cout << endl;
	}

	// ouput squared forces

	for (i=0; i<N; i++){		
		cout << "Squaed Force/Torque on particle i = " << i << ": ";
		cout << setw(30) << setprecision(10) << F[i][0]*F[i][0] + F[i][1]*F[i][1] + F[i][2]*F[i][2];	
		cout << setw(30) << setprecision(10) << TqW[i][0]*TqW[i][0] + TqW[i][1]*TqW[i][1] + TqW[i][2]*TqW[i][2];			
		cout << endl;
	}

	// get all perturbed forces
	for (j=0; j<N; j++){

		// loop over perturbation directions
		for (d2=0; d2<DOF; d2++){

			// get linear index
			l = DOF*j+d2;

			// perturb particle j in d2 (+h) direction, get updated forces
			this->rb_perturbation(j,d2,h);

			// loop over to extract forces to perturb
			for (i=0; i<N; i++){

				// perturb in DOF dimensions
				for (d1=0; d1<DOF; d1++){
					// get linear index
					k = DOF*i+d1;

					// add forces/torqu to matrix
					if (d1 < NDIM)
						Pplus(k,l) = F[i][d1];
					else
						Pplus(k,l) = TqW[i][d1];
				}
			}

			// perturb particle j in d2 (-h) direction, get updated forces
			this->rb_perturbation(j,d2,-h);

			// loop over to extract forces to perturb
			for (i=0; i<N; i++){

				// perturb in DOF dimensions
				for (d1=0; d1<DOF; d1++){
					// get linear index
					k = DOF*i+d1;

					// add forces/torques to matrix
					if (d1 < NDIM)
						Pminus(k,l) = F[i][d1];
					else
						Pminus(k,l) = TqW[i][d1];
				}
			}
		}
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
	obj << DOF << endl;
	obj << nr << endl;
	obj << ((N-nr)*DOF + NDIM - 1) - this->get_c_sum() << endl;
	for (k=0; k<DOF*N; k++){
		for (l=k; l<DOF*N; l++){
			// get matrix entry
			e = (N*DOF)*k + l - (((k+1)*k)/2);

			// calculate single matrix entry
			// Mkl = (Fk.at(k) - P(k,l) + Fk.at(l) - P(l,k))/(2*h);
			Mkl = (Pminus(k,l) - Pplus(k,l) + Pminus(l,k) - Pplus(l,k))/(4*h);
			dynMatrix(k,l) = Mkl;
			dynMatrix(l,k) = Mkl;

			// print matrix
			obj << setw(8) << k << setw(8) << l << setw(30) << setprecision(20) << Mkl << endl;
		}
	}

	// NOW print out mass matrix	
	for (i=0; i<N; i++){

		// first print out mass of residue
		obj << m[i] << endl;

		// counter
		k = DOF*i;

		// get first block in mass matrix
		for (d1=0; d1<NDIM; d1++){
			massMatrix(k+d1,k+d1) = m[i];
		}

		// now print out upper triangle of moment of inertia tensor IN WORLD FRAME
		for (d1=0; d1<NDIM; d1++){
			for (d2=d1; d2<NDIM; d2++){
				Iij = 0;
				if (d1==d2){
					// get diagonal element
					for (j=0; j<Na[i]; j++){
						// atomic mass
						am = (4.0/3.0) * PI * pow(ar[i][j],3);

						// lever arm length
						dnorm = sqrt(xW[i][j][0]*xW[i][j][0] + xW[i][j][1]*xW[i][j][1] + xW[i][j][2]*xW[i][j][2]);

						// MoI for sphere rotating around its own axis
						Sij = (2.0/5.0)*am*pow(ar[i][j],2);

						// Iij element from parallel axis theorem: https://en.wikipedia.org/wiki/Parallel_axis_theorem
						Iij += Sij + am*(dnorm - xW[i][j][d1]*xW[i][j][d1]);
					}
				}
				else{
					for (j=0; j<Na[i]; j++){
						// atomic mass
						am = (4.0/3.0) * PI * pow(ar[i][j],3);

						// Iij element from parallel axis theorem: https://en.wikipedia.org/wiki/Parallel_axis_theorem
						Iij -= am*xW[i][j][d1]*xW[i][j][d2];

					}
				}

				// print out Iij entry
				obj << Iij << endl;

				// save to massMatrix
				massMatrix(k+NDIM+d1,k+NDIM+d2) = Iij;
				massMatrix(k+NDIM+d2,k+NDIM+d1) = Iij;
			}
		}		
	}

	// do eigenvalue decomposition in c++
	cout << "doing eigenvalue decomposition" << endl;

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dynMatrix, massMatrix);
    cout << "Generalized eigenvalues (verify):\n" << es.eigenvalues() << endl;
    // cout << "Generalized eigenvectors (verify):\n" << es.eigenvectors().rightCols(N*DOF).topRows(N*DOF) << endl;

	obj.close();
}


// perturbation function
void rigidbody::rb_perturbation(int i, int d, double h){
	// local variables
	int i0,d0;
	double x0,theta0,dtheta;

	// reset force and torque
	for (i0=0; i0<N; i0++){
		for (d0=0; d0<NDIM; d0++){
			F[i0][d0] = 0.0;
			TqW[i0][d0] = 0.0;
		}

	}

	// perturb either translation or rotation
	if (d < NDIM){
		// translate particle i
		x0 = x[i][d];
		x[i][d] += h;

		// calc new force
		this->force_update();

		// reverse perturbation
		x[i][d] = x0;
	}
	else{		

		// get angle (path length of rotation/radius of rotater)
		dtheta = h/(0.77*r[i]);

		// rotate forward to test energy (d=3 => eulang1, d=4 => eulang2, etc)
		this->rotate_single_particle_xyz(i,d-3,dtheta);
		// this->rotate_single_particle(i,d-2,dtheta);

		// update forces based on new angle
		this->force_update();

		// rotate back to original location to reverse pertubation
		this->rotate_single_particle_xyz(i,d-3,-dtheta);
		// this->rotate_single_particle(i,d-2,-dtheta);
	}
}


// PE minimizer

void rigidbody::rb_fire_umin(int NTmax, double Ktol){
	// local variables
	int t,kr;	

	// constant energy checking
	double Uold, dU, dUtol;
	int epc, epcN, epconst;
	epconst = 0;
	Uold = 0;
	dU = 0;
	dUtol = 1e-12;
	epc = 0;
	epcN = 5e3;

	// initial K and U to 0
	K = 10*Ktol;
	U = 0;

	// start iterate at 0
	t = 0;
	
	// Minimize energy using FIRE
	while ((K > Ktol || epconst != 1) && t < NTmax) {
		// update nearest neighbor lists if applicable
		if (t % nnupdate == 0 && NCL > -1){
			cout << "^ ";
			this->update_nlcl(t);
		}

		// advance quaternions, positions
		this->verlet_first();

		// update forces
		this->force_update();

		// include fire relaxation
		this->rb_fire();

		// advance angular momentum
		this->verlet_second();	

		// check for rattlers
		kr = 0;
		nr = this->rmv_rattlers(kr);

		// check for constant potential energy
		dU = abs(Uold - U);
		if (dU < dUtol) {
			epc++;
			if (epc > epcN)
				epconst = 1;
		}
		else {
			epconst = 0;
			epc = 0;
		}
		Uold = U;	

		// output information
		if (t % plotskip == 0) {
			this->monitor_header(t);
			this->rigidbody_md_monitor();
			cout << "** ROOT SEARCH: " << endl;
			cout << "phi = " << phi << endl;
			cout << "epconst = " << epconst << endl;
			cout << "nr = " << nr << endl;
			cout << "alpha = " << alpha << endl;
			cout << endl;
			cout << endl;
		}

		// update iterate		
		t++;
	}
	if (t == NTmax){
		cout << "** COULD NOT FIND ENERGY MINIMUM IN NTmax TIME, ENDING WITH ERROR..." << endl;
		throw;
	}
}












void rigidbody::single_md(int t){
	// update nearest neighbor lists if applicable
	if (t % nnupdate == 0 && NCL > -1)
		this->update_nlcl(t);

	// advance quaternions, positions
	this->verlet_first();

	// update forces
	this->force_update();

	// advance angular momentum
	this->verlet_second();
}

void rigidbody::single_fire(int t){
	// update nearest neighbor lists if applicable
	if (t % nnupdate == 0 && NCL > -1)
		this->update_nlcl(t);

	// advance quaternions, positions
	this->verlet_first();

	// update forces
	this->force_update();

	// include fire relaxation
	this->rb_fire();

	// advance angular momentum
	this->verlet_second();
}





// FUNCTION to read in positions from config file,
//
//
//

// rigidbody::contact_break(int NT, double tmp0){

// }