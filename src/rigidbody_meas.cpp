/*

	Measure methods implementation
	for rigidbody class

	BY Jack Treado
	10/01/2018

*/

#include "rigidbody.h"
#include "vec3.h"
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


// functions for analytical dm
void rigidbody::populate_Pmatrix(Eigen::MatrixXd& P, int mu, int i){
	// populate matrix with position values
	P(0,0) = 0.0;
	P(1,0) = -xW[mu][i][2];
	P(2,0) = xW[mu][i][1];

	P(0,1) = xW[mu][i][2];
	P(1,1) = 0.0;
	P(2,1) = -xW[mu][i][0];

	P(0,2) = -xW[mu][i][1];
	P(1,2) = xW[mu][i][0];
	P(2,2) = 0.0;	
}

void rigidbody::calculate_vec_mat_cp(Eigen::MatrixXd& outMat, vec3& rightVec, vec3& matVec){
	// calculate column-wise cross product of vector with matrix (see math in Sec 4.5 of notes)

	// local variables
	int d1,d2;
	double vecDotP;

	// get dot product of vectors
	vecDotP = rightVec*matVec;

	// loop over matrix elements, compute outer product
	for (d1=0; d1<NDIM; d1++){
		for (d2=0; d2<NDIM; d2++){
			outMat(d1,d2) = rightVec.at(d1)*matVec.at(d2);
			if (d1 == d2)
				outMat(d1,d2) -= vecDotP;
		}
	}
}

void rigidbody::calculate_mass_matrix(Eigen::MatrixXd& massMatrix){
	// local variables
	int i,j,k,d1,d2;
	double Iij,Sij,am,dnorm;

	// NOW print out mass matrix	
	for (i=0; i<N; i++){

		// counter
		k = DOF*i;

		// get first block in mass matrix
		for (d1=0; d1<NDIM; d1++)
			massMatrix(k+d1,k+d1) = m[i];

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
					// get off-diagonal element
					for (j=0; j<Na[i]; j++){
						// atomic mass
						am = (4.0/3.0) * PI * pow(ar[i][j],3);

						// Iij element from parallel axis theorem: https://en.wikipedia.org/wiki/Parallel_axis_theorem
						Iij -= am*xW[i][j][d1]*xW[i][j][d2];

					}
				}

				// save to massMatrix
				massMatrix(k+NDIM+d1,k+NDIM+d2) = Iij;
				massMatrix(k+NDIM+d2,k+NDIM+d1) = Iij;
			}
		}		
	}
}

// Analytical form of the rigidbody dynamical matrix
void rigidbody::rb_analytical_dm(string& dmstr){
	// local variables
	int mu,nu,alpha,fd,qd;
	int i,j,k,l,kreal,lreal,k2real,l2real;
	double Rij,Sij,rij,sij,kij,hij;
	double aij[NDIM];

	// open dynamical matrix file
	ofstream obj(dmstr.c_str());
	if (!obj.is_open()){
		cout << "could not open file " << dmstr << endl;
		throw;
	}

	// matrices
	Eigen::MatrixXd M_translation_force(NDIM,NDIM);  		// matrix of derivatives of force wrt translations
	Eigen::MatrixXd M_rotation_force(NDIM,NDIM);			// matrix of derivatives of force wrt rotations
	Eigen::MatrixXd M_translation_torque(NDIM,NDIM);		// matrix of derivatives of torque wrt translations
	Eigen::MatrixXd M_rotation_torque(NDIM,NDIM);			// matrix of derivatives of torque wrt rotations
	Eigen::MatrixXd M_rotation_torque_diag(NDIM,NDIM);		// matrix of derivatives of torque wrt rotations: SPECIFICALLY FOR DIAGONAL CASE
	Eigen::MatrixXd P_mu_i(NDIM,NDIM);						// P matrix, made of cross product of r_mu_i with Id.
	Eigen::MatrixXd P_nu_j(NDIM,NDIM);						// P matrix, made of cross product of r_nu_j with Id.
	Eigen::MatrixXd cP_mu_i(NDIM,NDIM);						// matrix of rel pos of atom(mu,i) crossed into columns of P(mu,i) matrix
	Eigen::MatrixXd cP_nu_j(NDIM,NDIM);						// matrix of rel pos of atom(mu,i) crossed into columns of P(nu,j) matrix
	Eigen::MatrixXd cP_mu_ij(NDIM,NDIM);					// matrix of columns of P matrix crossed into rij connection vector 

	// dynamical matrix
	Eigen::MatrixXd dynMatrix(DOF*N,DOF*N);					// The Big Cheese Herself, the full dynamical matrix
	Eigen::MatrixXd massMatrix(DOF*N,DOF*N);				// The Big Cheese's little sister, the mass matrix

	// vec3 objects
	vec3 rij_uvec; 											// vector between atom (mu,i) and (nu,j)
	vec3 mu_rel_vec;										// position of (mu,i) relative to com
	vec3 nu_rel_vec;										// position of (nu,j) relative to com
	vec3 c_mu_ij;											// cross product with rel position of atom (mu,i) and rij
	vec3 c_nu_ij;											// cross product with rel position of atom (nu,j) and rij

	// loop over pairwise particles
	cout << "Looping over particles to update matrix entries..." << endl;
	for (mu=0; mu<N; mu++){

		// first loop over all particles, and calculate contribution 
		// to sum for diagonal block matrix
		for (nu=mu+1; nu<N; nu++){

			// check if in contact
			Rij = this->get_distance(mu,nu);

			// contact distance
			Sij = r[mu] + r[nu];

			// only check couplings if in contact
			if (Rij < Sij){

				// reset block matrices
				for (fd=0; fd<NDIM; fd++){
					for (qd=0; qd<NDIM; qd++){
						M_translation_force(fd,qd) = 0.0;
						M_rotation_force(fd,qd) = 0.0;
						M_translation_torque(fd,qd) = 0.0;
						M_rotation_torque(fd,qd) = 0.0;
						M_rotation_torque_diag(fd,qd) = 0.0;
					}
				}

				for (i=0; i<Na[mu]; i++){
					for (j=0; j<Na[nu]; j++){
						// distance between atoms
						rij = this->get_atomic_distance(mu,nu,i,j,aij);

						// contact distance
						sij = ar[mu][i] + ar[nu][j];						

						// only check coupling if atoms in contact
						if (rij < sij){

							// populate vec3 objects
							rij_uvec.reset(aij);
							rij_uvec.normalize();				
							mu_rel_vec.reset(xW[mu][i]);	
							nu_rel_vec.reset(xW[nu][i]);								

							// calculate cross products
							c_mu_ij = mu_rel_vec % rij_uvec;
							c_nu_ij = nu_rel_vec % rij_uvec;

							// populate P matrices		
							populate_Pmatrix(P_mu_i,mu,i);
							populate_Pmatrix(P_nu_j,nu,j);

							// calculate vector/matrix cross products
							calculate_vec_mat_cp(cP_nu_j,mu_rel_vec,nu_rel_vec);
							calculate_vec_mat_cp(cP_mu_i,mu_rel_vec,mu_rel_vec);
							calculate_vec_mat_cp(cP_mu_ij,rij_uvec,mu_rel_vec);							

							// calculate overlap
							hij = (1-(rij/sij));

							// effective spring constant
							kij = ep/(rij*sij);
							
							// fill block matrices (compute outer products)
							for (fd=0; fd<NDIM; fd++){
								for (qd=0; qd<NDIM; qd++){
									// 1. Translation-Force derivative
									M_translation_force(fd,qd) += kij*rij_uvec.at(fd)*rij_uvec.at(qd);
									if (fd == qd)
										M_translation_force(fd,qd) += kij*hij;

									// 2. Rotation-Force derivative
									M_rotation_force(fd,qd) += kij*(rij_uvec.at(fd)*c_nu_ij.at(qd) + hij*P_nu_j(fd,qd));

									// 3. Translation-Torque derivative
									M_translation_torque(fd,qd) += kij*(c_mu_ij.at(fd)*rij_uvec.at(qd) - hij*P_mu_i(fd,qd));

									// 4. Rotation-Torque derivatives (case where mu \neq nu)
									M_rotation_torque(fd,qd) += kij*(c_nu_ij.at(fd)*c_mu_ij.at(qd) + hij*cP_nu_j(fd,qd));

									// 4b. Rotation-Torque derivative (case where mu = nu)
									M_rotation_torque_diag(fd,qd) += kij*(c_mu_ij.at(fd)*c_mu_ij.at(qd) + hij*(cP_mu_i(fd,qd) - rij*cP_mu_ij(fd,qd)));
								}
							}
						}
					}
				}

				// update off-diagonal and diagonal blocks from sums over atoms

				// NOTE: since k,l both start at 0, this is a loop over the whole block matrix
				// D_mu_nu, so diagonal blocks are completely accounted for

				// ALSO NOTE: we access by k-NDIM and l-NDIM below when k,l > NDIM (for rotational dofs)
				// because above we allocate each block matrix as 3x3
				for (k=0; k<DOF; k++){
					for (l=0; l<DOF; l++){
						// entries in actual dynamical matrix
						kreal = DOF*mu + k;			// row of off-diagonal(upper T)
						lreal = DOF*nu + l;			// column of off-diagonal (upper T)
						k2real = DOF*mu + l;		// diagonal row
						l2real = DOF*nu + k;		// diagonal column

						// add elements from correct matrix, and calculate contribution to case when mu = nu
						if (k < NDIM){
							// translation-force
							if (l < NDIM){
								// off diagonal (= D_mu_nu)
								dynMatrix(kreal,lreal) = M_translation_force(k,l);

								// add to diagonal D_mu_mu
								dynMatrix(kreal,k2real) -= M_translation_force(k,l);

								// add TRANSPOSE to diagonal D_nu_nu
								dynMatrix(lreal,l2real) -= M_translation_force(k,l);
							}
							// rotation-force
							else{
								// off diagonal (= D_mu_nu)
								dynMatrix(kreal,lreal) = M_rotation_force(k,l-NDIM);

								// add to diagonal D_mu_mu
								dynMatrix(kreal,k2real) -= M_rotation_force(k,l-NDIM);

								// add TRANSPOSE of translation-torque to diagonal D_nu_nu
								dynMatrix(lreal,l2real) -= M_translation_torque(k,l-NDIM);
							}
						}
						else{
							// translation-torque
							if (l < NDIM){
								// off diagonal (= D_mu_nu)
								dynMatrix(kreal,lreal) = M_translation_torque(k-NDIM,l);

								// add to diagonal D_mu_mu
								dynMatrix(kreal,k2real) -= M_translation_torque(k-NDIM,l);

								// add TRANSPOSE of rotation-force to diagonal D_nu_nu
								dynMatrix(lreal,l2real) -= M_rotation_force(k-NDIM,l);
							}
							// rotation-torque
							else{
								// off diagonal (= D_mu_nu)
								dynMatrix(kreal,lreal) = M_rotation_torque(k-NDIM,l-NDIM);

								// add to diagonal D_mu_mu (slightly different than off-diag matrix)
								dynMatrix(kreal,k2real) -= M_rotation_torque_diag(k-NDIM,l-NDIM);

								// add TRANSPOSE to diagonal D_nu_nu (slightly different than off-diag matrix)
								dynMatrix(lreal,l2real) -= M_rotation_torque_diag(k-NDIM,l-NDIM);
							}
						}

						// make sure matrix is symmetric
						dynMatrix(lreal,kreal) = dynMatrix(kreal,lreal);
					}
				}				

			}
		}
	}

	// calculate mass matrix
	cout << "calculating mass matrix..." << endl;
	this->calculate_mass_matrix(massMatrix);

	// solve generalized eigenvalue equation
	cout << "doing eigenvalue decomposition" << endl;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dynMatrix, massMatrix);
    cout << "Generalized eigenvalues (verify):\n" << es.eigenvalues() << endl;

    // print eigenvalues of dynamical matrix info to file
    obj << N << endl;
	obj << DOF << endl;
	obj << nr << endl;
	obj << ((N-nr)*DOF + NDIM - 1) - this->get_c_sum() << endl;
	obj << dynMatrix << endl;
	obj << massMatrix << endl;
    obj << es.eigenvalues() << endl;

    // close object
    obj.close();
}


// numerical form
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
			Mkl = (Pminus(k,l) - Pplus(k,l) + Pminus(l,k) - Pplus(l,k))/(4.0*h);
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