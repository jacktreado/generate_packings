/*

	Methods implementation 
	for packing class

	BY Jack Treado

*/

#include "packing.h"

using namespace std;

const double PI = 4*atan(1);

// FUNCTION to populate entries in dynamical matrix
// assuming that forces are already calculated
void packing::dynamical_matrix(ofstream& dmobj, double h0){
	// local variables
	int i,d,nr,kr;
	double h,dt0,len;
	int Nentries = (NDIM*N*((NDIM*N)+1))/2; 			// number of unique entries in dynamical matrix	

	// local matrices
	Eigen::MatrixXd numHessian(NDIM*N,NDIM*N);			// numerical hessian
	Eigen::MatrixXd massMatrix(NDIM*N,NDIM*N);			// mass matrix
	Eigen::MatrixXd anHessian(NDIM*N,NDIM*N);			// analytical hessian

	// rescale all lengths to smaller radius
	len = r[0];
	ep = 1.0;
	this->rescale_lengths(len);

	// rescale time as well
	dt0 = 0.01;
	this->set_md_time(dt0);
	this->set_dtmax(10.0);

	// set h scale
	h = h0*L[0];

	// remove any rattlers
	kr = 0;
	nr = this->rmv_rattlers(kr);
	cout << "there are nr = " << nr << " rattlers" << endl;
	
	// print all unperturbed forces first
	cout << "getting forces..." << endl;	
	this->hs_force();
	cout << "total potential energy U = " << U << endl;
	for (i=0; i<N; i++){		
		cout << "Force on particle i = " << i << ": ";
		for (d=0; d<NDIM; d++)
			cout << setw(20) << setprecision(6) << F[i][d];
		cout << endl;
	}

	// compute numerical hessian and massMatrix
	this->compute_numerical_hessian(numHessian,massMatrix,h);

	// store all unperturbed forces first
	cout << "outputting forces again to check if perturbations remain..." << endl;
	U = 0;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			F[i][d] = 0;
	}
	this->hs_force();
	cout << "total potential energy U = " << U << endl;

	for (i=0; i<N; i++){		
		cout << "Force on particle i = " << i << ": ";
		for (d=0; d<NDIM; d++)
			cout << setw(20) << setprecision(6) << F[i][d];
		cout << endl;
	}	

	// print out matrix information
	dmobj << N << endl;
	dmobj << NDIM << endl;
	dmobj << nr << endl;

	// solve generalized eigenvalue problem for numerical hessian
	cout << "doing eigenvalue decomposition for numerical hessian" << endl;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> num_eigs(numHessian, massMatrix);
    cout << "Eigenvalues for numerical dynamical matrix (verify):\n" << num_eigs.eigenvalues() << endl;

    cout << "compare to analytical hessian:" << endl;
    compute_analytical_hessian(anHessian);
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> an_eigs(anHessian, massMatrix);
	cout << "Eigenvalues for analytical dynamical matrix:\n" << an_eigs.eigenvalues() << endl;

	// print everything to object
	dmobj << num_eigs.eigenvalues() << endl;
	dmobj << an_eigs.eigenvalues() << endl;
}

// compute analytical DM, excite all modes and measure VACF 
void packing::all_mode_perturbation(std::ofstream& dmobj, int NT, int vsave, double T0){
	// local variables
	int i,d,mtmp,NTeq,kr,nr;
	double wM,aM,pM,len,dt0,evecNorm,wLow;
	Eigen::VectorXd evecM(NDIM*N);

	// local matrices
	Eigen::MatrixXd massMatrix(NDIM*N,NDIM*N);			// mass matrix
	Eigen::MatrixXd hessian(NDIM*N,NDIM*N);				// analytical hessian

	// rescale all lengths to smaller radius
	len = r[0];
	ep = 1.0;
	this->rescale_lengths(len);

	// rescale time as well
	dt0 = 0.01;
	this->set_md_time(dt0);
	this->set_dtmax(10.0);

	// get rattlers
	kr = 0;
	nr = this->rmv_rattlers(kr);
	cout << "There are nr = " << nr << " rattlers, calculating analytical hessian " << endl;

	// compute values for eigenvalue decomposition
	compute_analytical_hessian(hessian);
	compute_mass_matrix(massMatrix);

	// solve generalized eigenvalue problem 
	cout << "Solving eigenvalues for analytical hessian" << endl;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigs(hessian, massMatrix);
	cout << "Generalized eigenvalues for analytical DM:\n" << eigs.eigenvalues() << endl;

	// find lowest non-trivial eigenvalue
	for (i=0; i<N*NDIM; i++){
		wLow = eigs.eigenvalues()[i];
		if (wLow > 1e-8){
			wLow = sqrt(wLow);
			cout << "smallest non-trivial eigenvalue = " << wLow << ", corresponding to mode = " << i << endl;
			break;
		}
	}

	// output initial info for dmobj
	dmobj << N << endl;
	dmobj << NDIM << endl;
	dmobj << nr << endl;
	dmobj << eigs.eigenvalues() << endl;

	// get NTeq (20 periods of lowest frequency eigenvalue)
	NTeq = 40*round((2*PI)/(wLow*dt));

	// calculate vacf
	cout << "calculating vacf with perturbation along all modes..." << endl;
	calc_vacf(NT,NTeq,vsave,1,T0,dmobj);
	cout << "vacf calculated, ending function single_mode_perturbation" << endl;
}

// compute analytical DM, excite mode and measure VACF to check for fidelity
void packing::single_mode_perturbation(std::ofstream& dmobj, int mode, int NT, int vsave, double T0){
	// local variables
	int i,d,mtmp,NTeq,kr,nr;
	double wM,aM,pM,len,dt0,evecNorm,wLow;
	Eigen::VectorXd evecM(NDIM*N);

	// local matrices
	Eigen::MatrixXd massMatrix(NDIM*N,NDIM*N);			// mass matrix
	Eigen::MatrixXd hessian(NDIM*N,NDIM*N);				// analytical hessian

	// rescale all lengths to smaller radius
	len = r[0];
	ep = 1.0;
	this->rescale_lengths(len);

	// rescale time as well
	dt0 = 0.01;
	this->set_md_time(dt0);
	this->set_dtmax(10.0);

	// get rattlers
	kr = 0;
	nr = this->rmv_rattlers(kr);
	cout << "There are nr = " << nr << " rattlers, calculating analytical hessian " << endl;

	// compute values for eigenvalue decomposition
	compute_analytical_hessian(hessian);
	compute_mass_matrix(massMatrix);

	// solve generalized eigenvalue problem 
	cout << "Solving eigenvalues for analytical hessian" << endl;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigs(hessian, massMatrix);
	cout << "Generalized eigenvalues for analytical DM:\n" << eigs.eigenvalues() << endl;

	// if mode -1, make sure wM is is lowest freq that is nontrivial
	mtmp = 0;
	for (i=0; i<NDIM*N; i++){
		// get frequency
		wM = abs(eigs.eigenvalues()[i]);

		if (wM < 1e-14)
			// skip trivial 0s
			continue;
		else{
			// save lowest frequency eigenmode
			if (mtmp == 0)
				wLow = sqrt(wM);

			// break when # of frequencies above trivial 0s = mode
			if (mtmp == mode)
				break;
			mtmp++;
		}
	}


	// select eigenvalue and eigenvector for mode
	mode = i;
	cout << "getting eigenvalue and eigenvector for mode = " << i << endl;
	wM = sqrt(eigs.eigenvalues()[mode]);
	evecM = eigs.eigenvectors().col(mode);

	// output initial info for dmobj
	dmobj << N << endl;
	dmobj << NDIM << endl;
	dmobj << nr << endl;
	dmobj << eigs.eigenvalues() << endl;
	dmobj << wM << endl;	

	// normalize evec
	evecNorm = 0.0;
	for (i=0; i<NDIM*N; i++)
		evecNorm += evecM(i)*evecM(i);

	evecNorm = sqrt(evecNorm);
	for (i=0; i<NDIM*N; i++)
		evecM(i) = evecM(i)/evecNorm;

	// initial conditions based on mode
	cout << "initial eigenvector" << endl;
	cout << evecM << endl;

	// eigenvector contribution to initial kinetic energy (pM)
	pM = 0.0;
	for (i=0; i<N; i++){

		// don't excite rattlers
		if (pc[i] == 0)
			continue;

		for (d=0; d<NDIM; d++)
			if (abs(evecM(NDIM*i+d)) > 1e-8)
				pM += evecM(NDIM*i+d)*evecM(NDIM*i+d)*m[i];
	}

	// amplitude to give kinetic energy T0
	aM = sqrt((2*T0)/pM);

	// set initial velocities
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			if (abs(evecM(NDIM*i+d)) < 1e-8)
				v[i][d] = 0.0;
			else if (pc[i] == 0)
				v[i][d] = 0.0;
			else
				v[i][d] = aM*evecM(NDIM*i+d);
		}
	}

	// get NTeq (20 periods of lowest frequency eigenvalue)
	NTeq = 40*round((2*PI)/(wLow*dt));

	// calculate vacf
	cout << "calculating vacf with perturbation along mode " << mode << endl;
	calc_vacf(NT,NTeq,vsave,0,T0,dmobj);
	cout << "vacf calculated, ending function single_mode_perturbation" << endl;
}

// numerical method to calculate the Hessian, to check mixed numerical method
void packing::compute_numerical_hessian(Eigen::MatrixXd& hessian, Eigen::MatrixXd& massMatrix, double h){
	// local variables
	int i,j,k,l,d1,d2;
	double Mkl;

	// local matrices
	Eigen::MatrixXd Pplus(NDIM*N,NDIM*N);				// matrix of all perturbed force (+h)
	Eigen::MatrixXd Pminus(NDIM*N,NDIM*N);				// matrix of all perturbed force (-h)

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

	// loop over particles, get forces as a function of particle perturbations
	// definition from: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
	for (k=0; k<NDIM*N; k++){
		for (l=k; l<NDIM*N; l++){
			// calculate single matrix entry
			Mkl = (Pminus(k,l) - Pplus(k,l) + Pminus(l,k) - Pplus(l,k))/((4.0/5.0)*h);

			// store matrix element in dynamical matrix
			hessian(k,l) = Mkl;
			hessian(l,k) = Mkl;
		}

		// also populate mass matrix
		i = k % NDIM;
		massMatrix(k,k) = m[i];
	}
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

// analytical method to calculate the Hessian, to check mixed numerical method
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

		if (pc[i] == 0){
			for (l=k; l<DF; l++){
				hessian(k,l) = 0.0;
				hessian(l,k) = 0.0;
			}
			continue;
		}

		// get dimension d1
		d1 = k % NDIM;
		for (l=k; l<DF; l++){
			// get particle j
			j = floor(l/NDIM);

			if (pc[j] == 0){
				hessian(k,l) = 0.0;
				hessian(l,k) = 0.0;
				continue;
			}

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

// compute mass matrix real quick
void packing::compute_mass_matrix(Eigen::MatrixXd& massMatrix){
	int k,i;

	for (k=0; k<NDIM*N; k++){
		// get particle index
		i = k % NDIM;

		// fill matrix
		massMatrix(k,k) = m[i];
	}
}




// calculate vacf from already initialized system
void packing::calc_vacf(int NT, int NTeq, int vsave, int vinit, double T0, ofstream& obj){
	// get number of time steps
	int i,d,t,tcurrent,nr,kr;

	if (vinit == 1)
		this->rand_vel_init(T0);
	else
		cout << "running vacf with pre-set velocities " << endl;

	// get rattlers
	kr = 0;
	nr = this->rmv_rattlers(kr);
	
	// initialize vacf arrays
	int nsamp = NT/vsave;
	vector<double> vacf;
	vector<double> numer;
	double denom = 0.0;
	for (i=0; i<nsamp; i++){
		numer.push_back(0.0);
		vacf.push_back(0.0);
	}

	// initialize vacf list (vlist)
	vlist = new vector<double>*[N];
	for (i=0; i<N; i++){
		// allocate memory for single particle
		vlist[i] = new vector<double>[NDIM];

		// setup each vector
		for (d=0; d<NDIM; d++)
			vlist[i][d].resize(nsamp,0.0);
	}

	// initialize saving time index
	tcurrent = 0;

	// eqilibrate
	cout << "running equilibration before VACF data sampling";
	for (t=0; t<NTeq; t++){
		// first md step
		this->pos_update();

		// force update
		this->hs_force();

		// vel update (neglect rattlers)
		this->vel_update(1);

		// print info
		if (t % 10000 == 0){
			cout << ".";
		}
	}

	// loop over time
	cout << endl << endl << "looping over time, saving velocities" << endl;
	for (t=0; t<NT; t++){
		// first md step
		this->pos_update();

		// force update
		this->hs_force();

		// vel update (neglect rattlers)
		this->vel_update(1);

		// print info
		if (t % plotskip == 0){
			this->md_monitor(t,nr,-1,-1);
			cout << endl << endl;
			if (xyzobj.is_open())
				this->print_xyz();
		}

		// get vacf at right time
		if (t % vsave == 0){
			this->get_vacf(tcurrent,numer,denom);
			tcurrent++;
		}
	}

	// finish off the vacf!
	this->finish_vacf(NT,vsave,numer,denom,vacf);

	cout << "VACF MD completed!" << endl;
	this->print_vacf(NT,vsave,vacf,obj);
}

// get vacf from previous velocities
void packing::get_vacf(int tcurrent, vector<double>& numer, double& denom){	
	int i,j,d,t0,sampi;
	double v0,vplusdt;

	// push back current velocities
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			vlist[i][d].at(tcurrent) = v[i][d];
	}

	// loop over previous velocities, calc vacf
	for (j=0; j<=tcurrent; j++){
		// get t0, time between two points
		t0 = tcurrent - j;

		// add to numerator
		for (i=0; i<N; i++){
			// skip rattlers
			if (pc[i] == 0)
				continue;

			// loop over velocities, add to VACF
			for (d=0; d<NDIM; d++){
				// velocity at current time t
				v0 = vlist[i][d].at(t0);

				// velocity at time t + dt
				vplusdt = vlist[i][d].at(tcurrent);

				// add to numerator
				numer.at(j) += v0*vplusdt;
			}
		}
	}

	// add to denominator
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			denom += v[i][d]*v[i][d];
	}	
}

void packing::finish_vacf(int NT, int vsave, vector<double>& numer, double& denom, vector<double>& vacf){
	int i,j,d,nsamp;

	// // average denominator
	// denom /= NT*N;

	// loop over previous velocities, calc vacf
	nsamp = NT/vsave;

	for (j=0; j<nsamp; j++){
		// // average numerator term
		// numer.at(j) /= NT*N;

		// get quotient
		vacf.at(j) = numer.at(j)/denom;			
	}
}

void packing::print_vacf(int NT, int vsave, vector<double>& vacf, ofstream& obj){
	int s = vacf.size();
	int i,nsamp;
	nsamp = NT/vsave;

	obj << NT << endl;
	obj << nsamp << endl;
	obj << vsave << endl;
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








