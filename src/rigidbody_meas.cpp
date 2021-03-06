/*

	Measure methods implementation
	for rigidbody class

	BY Jack Treado
	10/01/2018

*/

#include "rigidbody.h"

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
void rigidbody::compute_rb_analytical_hessian(Eigen::MatrixXd& hessian){
	// local variables
	int mu,nu,alpha,fd,qd;
	int i,j,k,l,kreal,lreal,k2real,l2real;
	double Rij,Sij,rij,sij,kij,hij;
	double aij[NDIM];

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
			Rij = get_distance(mu,nu);

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
						rij = get_atomic_distance(mu,nu,i,j,aij);

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
								hessian(kreal,lreal) = M_translation_force(k,l);

								// add to diagonal D_mu_mu
								hessian(kreal,k2real) -= M_translation_force(k,l);

								// add TRANSPOSE to diagonal D_nu_nu
								hessian(lreal,l2real) -= M_translation_force(k,l);
							}
							// rotation-force
							else{
								// off diagonal (= D_mu_nu)
								hessian(kreal,lreal) = M_rotation_force(k,l-NDIM);

								// add to diagonal D_mu_mu
								hessian(kreal,k2real) -= M_rotation_force(k,l-NDIM);

								// add TRANSPOSE of translation-torque to diagonal D_nu_nu
								hessian(lreal,l2real) -= M_translation_torque(k,l-NDIM);
							}
						}
						else{
							// translation-torque
							if (l < NDIM){
								// off diagonal (= D_mu_nu)
								hessian(kreal,lreal) = M_translation_torque(k-NDIM,l);

								// add to diagonal D_mu_mu
								hessian(kreal,k2real) -= M_translation_torque(k-NDIM,l);

								// add TRANSPOSE of rotation-force to diagonal D_nu_nu
								hessian(lreal,l2real) -= M_rotation_force(k-NDIM,l);
							}
							// rotation-torque
							else{
								// off diagonal (= D_mu_nu)
								hessian(kreal,lreal) = M_rotation_torque(k-NDIM,l-NDIM);

								// add to diagonal D_mu_mu (slightly different than off-diag matrix)
								hessian(kreal,k2real) -= M_rotation_torque_diag(k-NDIM,l-NDIM);

								// add TRANSPOSE to diagonal D_nu_nu (slightly different than off-diag matrix)
								hessian(lreal,l2real) -= M_rotation_torque_diag(k-NDIM,l-NDIM);
							}
						}

						// make sure matrix is symmetric
						hessian(lreal,kreal) = hessian(kreal,lreal);
					}
				}				

			}
		}
	}	
}


// numerical form
void rigidbody::rb_dynamical_matrix(string& dmstr, double h0){
	// local variables
	int i,j,d,k,l,e,d1,d2,nr,kr;
	double Mkl,dMk,dMl,h,Sij,Iij,am,dnorm,rmin;	

	// rescale system to H atom = 1.0 A, unit energy scale
	rmin = 1e6;
	for (i=0; i<Na[0]; i++){
		if (ar[0][i] < rmin)
			rmin = ar[0][i];
	}
	rb_rescale(rmin);
	set_ep(1.0);
	h = h0*L[0];

	// remove any rattlers
	kr = 0; 
	nr = rmv_rattlers(kr);
	cout << "there are nr = " << nr << " rattlers" << endl;

	int Nentries = (DOF*N*((DOF*N)+1))/2; 			// number of unique entries in dynamical matrix
	Eigen::MatrixXd Pplus(DOF*N,DOF*N);				// matrix of all perturbed force (+h)
	Eigen::MatrixXd Pminus(DOF*N,DOF*N);			// matrix of all perturbed force (-h)

	// matrices to do eigenvalue decomposition in c++
	Eigen::MatrixXd numHessian(DOF*N,DOF*N);
	Eigen::MatrixXd anHessian(DOF*N,DOF*N);
	Eigen::MatrixXd massMatrix(DOF*N,DOF*N);	

	// output all forces to check
	cout << "outputting forces..." << endl;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			F[i][d] = 0.0;
			TqW[i][d] = 0.0;
		}
	}
	force_update();
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
		cout << "Magnitude of Force/Torque on particle i = " << i << ": ";
		cout << setw(30) << setprecision(10) << sqrt(F[i][0]*F[i][0] + F[i][1]*F[i][1] + F[i][2]*F[i][2]);	
		cout << setw(30) << setprecision(10) << sqrt(TqW[i][0]*TqW[i][0] + TqW[i][1]*TqW[i][1] + TqW[i][2]*TqW[i][2]);			
		cout << endl;
	}

	// calculate arc rad, i.e. mean residue length scale
	double arcrad = 0.0;
	for (i=0; i<N; i++)
		arcrad += r[i];
	arcrad /= N;

	// ONLY PERTURB FORCES THAT ARE PART OF THE ORIGINAL CONTACT MATRIX

	// get all perturbed forces
	for (j=0; j<N; j++){

		// loop over perturbation directions
		for (d2=0; d2<DOF; d2++){

			// get linear index
			l = DOF*j+d2;

			// perturb particle j in d2 (+h) direction, get updated forces
			rb_perturbation(j,d2,h,arcrad);

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
						Pplus(k,l) = TqW[i][d1-NDIM];
				}
			}

			if (d2 == 0 || d2 == 3){
				cout << endl << "New force and torque for j pert = "  << j << ", and d2 = " << d2 << endl;
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
			}

			// perturb particle j in d2 (-h) direction, get updated forces
			rb_perturbation(j,d2,-h,arcrad);

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
						Pminus(k,l) = TqW[i][d1-NDIM];
				}
			}
		}
	}


	// loop over particles, get forces as a function of particle perturbations
	// definition from: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm	

	// calculate numerical matrix elements
	// for (i=0; i<N; i++){
	// 	for (j=i; j<N; j++){
	// 		for (d1=0; d1<DOF; d1++){
	// 			for (d2=d1; d2<DOF; d2++){
	// 				// get matrix indices
	// 				k = DOF*i + d1;
	// 				l = DOF*j + d2;

	// 				// calculate different deriatives depending upon whether or not force or torque
	// 				if (d1 < NDIM)
	// 					dMk = (Pminus(k,l) - Pplus(k,l))/(4.0*h);
	// 				else
	// 					dMk = (Pminus(k,l) - Pplus(k,l))/(4.0*h);

	// 				if (d2 < NDIM)
	// 					dMl = (Pminus(l,k) - Pplus(l,k))/(4.0*h);
	// 				else
	// 					dMl = (Pminus(l,k) - Pplus(l,k))/(4.0*h);

	// 				// calc matrix element
	// 				Mkl = dMk + dMl;

	// 				// store in matrix, enfore symmetry
	// 				numHessian(k,l) = Mkl;
	// 				numHessian(l,k) = Mkl;
	// 			}
	// 		}
	// 	}
	// }

	for (k=0; k<DOF*N; k++){
		for (l=k; l<DOF*N; l++){
			// calculate single matrix entry
			// Mkl = (Fk.at(k) - P(k,l) + Fk.at(l) - P(l,k))/(2*h);
			Mkl = (Pminus(k,l) - Pplus(k,l) + Pminus(l,k) - Pplus(l,k))/(4.0*h);
			numHessian(k,l) = Mkl;
			numHessian(l,k) = Mkl;
		}
	}

	// calculate mass matrix
	cout << "calculating mass matrix..." << endl;
	calculate_mass_matrix(massMatrix);

	// solve generalized eigenvalue equation
	cout << "doing eigenvalue decomposition on numerical Hessian" << endl;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> numEigs(numHessian, massMatrix);
    cout << "Generalized numerical eigenvalues (verify):\n" << numEigs.eigenvalues() << endl;

    // // get analytical hessian    
    // cout << "computing elements of analytical hessian matrix..." << endl;
    // compute_rb_analytical_hessian(anHessian);

    // cout << "doing eigenvalue decomposition on analytical Hessian" << endl;
    // Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> anEigs(anHessian, massMatrix);

    // print dynamical matrix to file
	ofstream obj(dmstr.c_str());
	if (!obj.is_open()){
		cout << "could not open file " << dmstr << endl;
		throw;
	}

    // print eigenvalues of dynamical matrix info to file
    obj << N << endl;
	obj << DOF << endl;
	obj << nr << endl;
    obj << numEigs.eigenvalues() << endl;

    // loop over modes, print eigenvectors
    for (k=0; k<DOF*N; k++)
    	obj << numEigs.eigenvectors().col(k) << endl;

    // close object
    obj.close();
}

void rigidbody::dimer_dynamical_matrix(string& dmstr, double h0){
	// local variables
	int i,j,d,k,l,e,d1,d2,nr,kr;
	double Mkl,dMk,dMl,h,Sij,Iij,am,dnorm,rmin;	

	// rescale system to H atom = 1.0 A, unit energy scale
	rmin = 1e6;
	for (i=0; i<Na[0]; i++){
		if (ar[0][i] < rmin)
			rmin = ar[0][i];
	}
	rb_rescale(rmin);
	set_ep(1.0);
	h = h0*L[0];

	// remove any rattlers
	kr = 0; 
	nr = rmv_rattlers(kr);
	cout << "there are nr = " << nr << " rattlers" << endl;

	int Nentries = (DOF*N*((DOF*N)+1))/2; 			// number of unique entries in dynamical matrix
	Eigen::MatrixXd Pplus(DOF*N,DOF*N);				// matrix of all perturbed force (+h)
	Eigen::MatrixXd Pminus(DOF*N,DOF*N);			// matrix of all perturbed force (-h)

	// matrices to do eigenvalue decomposition in c++
	Eigen::MatrixXd numHessian(DOF*N,DOF*N);
	Eigen::MatrixXd massMatrix(DOF*N,DOF*N);

	// output all forces to check
	cout << "outputting forces..." << endl;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			F[i][d] = 0.0;
			TqW[i][d] = 0.0;
		}
	}
	force_update();
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
		cout << "Magnitude of Force/Torque on particle i = " << i << ": ";
		cout << setw(30) << setprecision(10) << sqrt(F[i][0]*F[i][0] + F[i][1]*F[i][1] + F[i][2]*F[i][2]);	
		cout << setw(30) << setprecision(10) << sqrt(TqW[i][0]*TqW[i][0] + TqW[i][1]*TqW[i][1] + TqW[i][2]*TqW[i][2]);			
		cout << endl;
	}

	// calculate arc rad, i.e. mean residue length scale
	double arcrad = 0.0;
	for (i=0; i<N; i++)
		arcrad += r[i];
	arcrad /= N;

	// ONLY PERTURB FORCES THAT ARE PART OF THE ORIGINAL CONTACT MATRIX

	// get all perturbed forces
	for (j=0; j<N; j++){

		// loop over perturbation directions
		for (d2=0; d2<DOF; d2++){

			// get linear index
			l = DOF*j+d2;

			// perturb particle j in d2 (+h) direction, get updated forces
			dimer_perturbation(j,d2,h,arcrad);

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
						Pplus(k,l) = TqW[i][d1-NDIM];
				}
			}

			if (d2 == 0 || d2 == 3){
				cout << endl << "New force and torque for j pert = "  << j << ", and d2 = " << d2 << endl;
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
			}

			// perturb particle j in d2 (-h) direction, get updated forces
			dimer_perturbation(j,d2,-h,arcrad);

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
						Pminus(k,l) = TqW[i][d1-NDIM];
				}
			}
		}
	}


	// loop over particles, get forces as a function of particle perturbations
	// definition from: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm	

	// calculate numerical matrix elements
	for (i=0; i<N; i++){
		for (j=i; j<N; j++){
			for (d1=0; d1<DOF; d1++){
				for (d2=d1; d2<DOF; d2++){
					// get matrix indices
					k = DOF*i + d1;
					l = DOF*j + d2;

					// calculate different deriatives depending upon whether or not force or torque
					if (d1 < NDIM)
						dMk = (Pminus(k,l) - Pplus(k,l))/(4.0*h);
					else
						dMk = (Pminus(k,l) - Pplus(k,l))/(4.0*h);

					if (d2 < NDIM)
						dMl = (Pminus(l,k) - Pplus(l,k))/(4.0*h);
					else
						dMl = (Pminus(l,k) - Pplus(l,k))/(4.0*h);

					// calc matrix element
					Mkl = dMk + dMl;

					// store in matrix, enfore symmetry
					numHessian(k,l) = Mkl;
					numHessian(l,k) = Mkl;

					// also enter values into mass matrix
					if (d1 == d2){
						if (d1 < NDIM)
							massMatrix(k,l) = m[i];
						else
							massMatrix(k,l) = Inn[i][d1 - NDIM];
					}
					else{
						massMatrix(k,l) = 0.0;
						massMatrix(l,k) = 0.0;
					}
				}
			}
		}
	}

	// solve generalized eigenvalue equation
	cout << "doing eigenvalue decomposition on numerical Hessian" << endl;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> numEigs(numHessian, massMatrix);
    cout << "Generalized numerical eigenvalues (verify):\n" << numEigs.eigenvalues() << endl;

    // print dynamical matrix to file
	ofstream obj(dmstr.c_str());
	if (!obj.is_open()){
		cout << "could not open file " << dmstr << endl;
		exit(1);
	}

    // print eigenvalues of dynamical matrix info to file
    obj << N << endl;
	obj << DOF << endl;
	obj << nr << endl;
    obj << numEigs.eigenvalues() << endl;

    // close object
    obj.close();
}


// perturbation functions
void rigidbody::rb_perturbation(int i, int d, double h, double arcrad){
	// local variables
	int i0,a0,d0;
	double x0,theta0,dtheta,radtmp;

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

		// calc new force ONLY FOR CURRENT CONTACT MATRIX
		contact_forces();

		// reverse perturbation
		x[i][d] = x0;
	}
	else{
		// store initial positions
		vector< vector<double> > x0;
		vector<double> xtmp(NDIM,0.0);
		for (a0=0; a0<Na[i]; a0++){
			for (d0=0; d0<NDIM; d0++)
				xtmp.at(d0) = xW[i][a0][d0];
			x0.push_back(xtmp);
		}

		// get angle (path length of rotation/radius of rotater)
		// dtheta = h/arcrad;
		dtheta = h;

		// rotate forward to test energy (d=3 => eulang1, d=4 => eulang2, etc)
		rotate_single_particle_xyz(i,d-3,dtheta);
		// rotate_single_particle(i,d-2,dtheta);

		// calc new force ONLY FOR CURRENT CONTACT MATRIX
		contact_forces();

		// // rotate back to original location to reverse pertubation
		// rotate_single_particle_xyz(i,d-3,-dtheta);
		for (a0=0; a0<Na[i]; a0++){
			for (d0=0; d0<NDIM; d0++)
				xW[i][a0][d0] = x0.at(a0).at(d0);
		}
	}
}

void rigidbody::dimer_perturbation(int i, int d, double h, double arcrad){
	// local variables
	int i0,a0,d0;
	double x0,theta0,dtheta,radtmp;

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

		// calc new force ONLY FOR CURRENT CONTACT MATRIX
		contact_forces();

		// reverse perturbation
		x[i][d] = x0;
	}
	else{
		// store initial positions
		vector< vector<double> > x0;
		vector<double> xtmp(NDIM,0.0);
		for (a0=0; a0<Na[i]; a0++){
			for (d0=0; d0<NDIM; d0++)
				xtmp.at(d0) = xW[i][a0][d0];
			x0.push_back(xtmp);
		}

		// get angle (path length of rotation/radius of rotater)
		dtheta = h/arcrad;

		// rotate forward to test energy
		rigidbody_xyz();
		if (d == NDIM)
			rotate_single_particle_xyz(i,1,dtheta);
		else if (d == NDIM + 1)
			rotate_single_particle_xyz(i,2,dtheta);
		else{
			cout << "	** ERROR: in dimer_perturbation(), cannot pass d=" << d << " to dimer pertubation, can only be NDIM or NDIM+1. Ending" << endl;
			exit(1);
		}
		rigidbody_xyz();

		// calc new force ONLY FOR CURRENT CONTACT MATRIX
		contact_forces();

		// // rotate back to original location to reverse pertubation
		// for (a0=0; a0<Na[i]; a0++){
		// 	for (d0=0; d0<NDIM; d0++)
		// 		xW[i][a0][d0] = x0.at(a0).at(d0);
		// }
		// // rigidbody_xyz();
		rotate_single_particle(i,d-2,-dtheta);
		rigidbody_xyz();
	}
}


// update forces solely based on contact network
void rigidbody::contact_forces(){
	// local variables
	int i, j, cj, ai, aj, d;
	double sij, rij, dx, da;
	double qix, qiy, qiz, qjx, qjy, qjz;
	double tix, tiy, tiz, tjx, tjy, tjz;
	double aij[NDIM];
	double fij[NDIM];

	// loop over contact matrix, calculate all forces based on contact network
	for (i=0; i<N; i++){

		// skip particles with no contacts
		if (pc[i] == 0)
			continue;

		for (j=i+1; j<N; j++){

			// reset torque counters
			tix = 0; tiy = 0; tiz = 0;
			tjx = 0; tjy = 0; tjz = 0;

			// skip particles with no contacts
			if (pc[i] == 0)
				continue;

			// mapping from matrix space to sub matrix space
			cj = N * i + j - ((i + 1) * (i + 2)) / 2;	

			// if in contact, calculate forces
			if (c[cj] == 1){
				for (ai = 0; ai < Na[i]; ai++) {
					for (aj = 0; aj < Na[j]; aj++) {
						// contact distance
						rij = ar[i][ai] + ar[j][aj];

						// get distance between atoms
						da = get_atomic_distance(i, j, ai, aj, aij);

						// if true, atoms are overlapping, so calc force and torquez
						if (da < rij) {
							// update contact forces
							for (d = 0; d < NDIM; d++) {
								fij[d] = hs(rij, da) * aij[d];
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




// PE minimizer
void rigidbody::rb_fire_umin(int NTmax, double Ktol){
	// local variables
	int i,d,t,kr,acsum,pcsum;	

	// constant energy checking
	double Uold, dU, dUtol, Fmag;
	int epc, epcN, epconst;
	epconst = 0;
	Uold = 0;
	dU = 0;
	dUtol = 1e-8;
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
			update_nlcl(t);
		}

		// advance quaternions, positions
		verlet_first();

		// update forces
		force_update();

		// include fire relaxation
		rb_fire();

		// advance angular momentum
		verlet_second();	

		// check for rattlers
		kr = 0;
		nr = rmv_rattlers(kr);

		// check for constant potential energy
		dU = abs((Uold - U)/Uold);
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
			monitor_header(t);
			rigidbody_md_monitor();
			cout << "** ROOT SEARCH: " << endl;
			cout << "phi = " << phi << endl;
			cout << "epconst = " << epconst << endl;
			cout << "nr = " << nr << endl;
			cout << "alpha = " << alpha << endl;
			cout << "dU = " << dU << endl;

			Fmag = 0.0;
			for (i=0; i<N; i++){
				for (d=0; d<NDIM; d++){
					Fmag += F[i][d]*F[i][d];
					Fmag += TqW[i][d]*TqW[i][d];
				}
			}
			Fmag = sqrt(Fmag);
			cout << "Fmag = " << Fmag << endl;

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
	else{
		acsum = 0.5*get_ac_sum();
		pcsum = get_c_sum();

		cout << "Energy minization completed!" << endl;
		cout << "Final U = " << U << endl;
		cout << "Final K = " << K << endl;
		cout << "Final phi = " << setprecision(6) << phi << endl;
		cout << "Final pcsum = " << pcsum << endl;
		cout << "Final acsum = " << acsum << endl;
		cout << "Final niso = " << (N-nr)*DOF - NDIM + 1 << endl;
		cout << "Final niso max = " << DOF*N - NDIM + 1 << endl;
		cout << "Final rattler # = " << nr << endl;
		cout << "Final contacts:" << endl;
		cout << "pc:" << endl;
		for (int i = 1; i < N + 1; i++) {
			cout << setw(6) << pc[i - 1];
			if (i % 10 == 0)
				cout << endl;
		}
		cout << endl;
		cout << "ac:" << endl;
		for (int i = 1; i < N + 1; i++) {
			cout << setw(6) << ac[i - 1];
			if (i % 10 == 0)
				cout << endl;
		}
		cout << endl;
		if (N < 40) {
			cout << "Contact Matrix:" << endl;
			print_c_mat();
		}
		if (xyzobj.is_open())
			rigidbody_xyz();
	}
}



// out velocities for res VACF trajectory
void rigidbody::rb_md_velocity(double T0, int NT, int tsave, ofstream& obj){
	// local variables
	int i,t,NTeq,kr,nr;
	double rmin;
	bool neglect_rattlers = true;

	// rescale system to H atom = 1.0 A, unit energy scale
	ep = 1.0;
	rmin = 1e6;
	for (i=0; i<Na[0]; i++){
		if (ar[0][i] < rmin)
			rmin = ar[0][i];
	}
	rb_rescale(rmin);
	set_md_time(0.01);
	set_dtmax(10.0);


	// remove any rattlers
	kr = 0; 
	nr = rmv_rattlers(kr);
	cout << "Before velocity run: There are nr = " << nr << " rattlers" << endl;

	// randomly assign velocities
	rand_vel_init(T0);

	// output header values to object
	obj << N << endl;			// number of residues
	obj << T0 << endl;			// constant temperature
	obj << NT << endl;			// number of time steps	
	obj << tsave << endl;		// time between samples
	obj << dt << endl;			// MD time scale in particle units

	// loop over over equilibration time (always 1 quarter of NT)
	NTeq = NT/4;

	for (t=0; t<NTeq; t++){
		// update nearest neighbor lists if applicable
		if (t % nnupdate == 0 && NCL > -1){
			cout << "^ ";
			update_nlcl(t);
		}

		// advance quaternions, positions
		verlet_first();

		// update forces
		force_update();

		// advance angular momentum
		verlet_second(neglect_rattlers);

		if (t % (5*plotskip) == 0){
			cout << " EQUILIBRATION " << endl;
			monitor_header(t);
			rigidbody_md_monitor();			
			cout << endl;
			cout << endl;
		}

		// rescale velocities to maintain temperature
		rescale_velocities(T0);
	}

	// loop over time, output variables every vsave
	for (t=0; t<NT; t++){
		// update nearest neighbor lists if applicable
		if (t % nnupdate == 0 && NCL > -1){
			cout << "^ ";
			update_nlcl(t);
		}

		// advance quaternions, positions
		verlet_first();

		// update forces
		force_update();

		// advance angular momentum
		verlet_second(neglect_rattlers);

		// output trajectory every vsave steps
		if (t % tsave == 0)
			output_velocities(obj);

		if (t % plotskip == 0){
			cout << " PRODUCTION " << endl;
			monitor_header(t);
			rigidbody_md_monitor();
			cout << endl;
			cout << endl;
		}
	}
}

// output velocitie and angular velocities for every particle
void rigidbody::output_velocities(ofstream& obj){
	// local variables
	int i,d,w,p;

	// width and precision of output
	w = 35;
	p = 16;

	// output particle velocities and angular velocities
	for (i=0; i<N; i++){
		// velocities
		for (d=0; d<NDIM; d++)
			obj << setw(w) << setprecision(p) << v[i][d];		

		// angular velocities
		for (d=0; d<NDIM; d++)
			obj << setw(w) << setprecision(p) << wM[i][d];

		// line break for next particle
		obj << endl;
	}
}



void rigidbody::single_md(int t){
	// update nearest neighbor lists if applicable
	if (t % nnupdate == 0 && NCL > -1)
		update_nlcl(t);

	// advance quaternions, positions
	verlet_first();

	// update forces
	force_update();

	// advance angular momentum
	verlet_second();
}

void rigidbody::single_fire(int t){
	// update nearest neighbor lists if applicable
	if (t % nnupdate == 0 && NCL > -1)
		update_nlcl(t);

	// advance quaternions, positions
	verlet_first();

	// update forces
	force_update();

	// include fire relaxation
	rb_fire();

	// advance angular momentum
	verlet_second();
}

