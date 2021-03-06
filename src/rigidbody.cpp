/*

	Methods implementation
	for rigidbody class

	BY Jack Treado

*/

#include "rigidbody.h"

using namespace std;

const double PI = atan(1)*4;

/*
==================================

	CONSTRUCTORS & DESTRUCTORS

==================================
*/



// constructor for residues to be read in
rigidbody::rigidbody(string &rbstr, int n, int dof, int nc, int s) : packing(n, dof, nc, s) {
	cout << "Entering rigidbody constructor..." << endl;
	cout << endl;

	cout << "N = " << N << endl;
	cout << "dof = " << dof << endl;
	cout << "nc = " << nc << endl;
	cout << "seed = " << seed << endl;

	// read in particle information
	cout << "getting info from file header" << endl;
	get_file_header(rbstr);

	// initialize particle info
	cout << "initializing particle variables" << endl;
	initialize_particles();

	// initialize dynamics
	cout << "initializing dynamics" << endl;
	initialize_dynamics();

	// read in rigid body info
	cout << "reading in data from " << rbstr << endl;
	read_in_info(rbstr);

	// get initial packing fraction	
	update_phi();
	cout << "Original phi = " << phi << endl;

	// initialize quaternions
	cout << "initializing quaternions" << endl;
	initialize_quaternions();

	// setup neighbor & cell list if nc is positive
	if (nc > -1){
		cout << "initializing nlcl" << endl;
		initialize_nlcl();
	}
	else
		cout << "nc = -1, no nlcl for this simulation" << endl;

	cout << endl;
	cout << "... rigidbody constructor complete." << endl;
	cout << endl << endl << endl;
}

rigidbody::~rigidbody() {
	cout << "~entering rigidbody destructor" << endl;

	// delete double & triple arrays
	int i, j;

	for (i = 0; i < N; i++) {
		// delete extra dimensions
		for (j = 0; j < Na[i]; j++) {
			delete [] xW[i][j];
			delete [] xM[i][j];
		}

		// triple arrays
		delete [] xW[i];
		delete [] xM[i];

		// double arrays
		delete [] wW[i];
		delete [] wM[i];

		delete [] LW[i];
		delete [] LWhalf[i];
		delete [] LM[i];
		delete [] LMhalf[i];
		delete [] LMdot[i];

		delete [] TqW[i];
		delete [] TqM[i];

		delete [] Inn[i];

		delete [] ar[i];
	}

	// delete outer shell pointers
	delete [] xW;
	delete [] xM;
	delete [] wW;
	delete [] wM;
	delete [] LW;
	delete [] LWhalf;
	delete [] LM;
	delete [] LMhalf;
	delete [] LMdot;
	delete [] TqW;
	delete [] TqM;
	delete [] Inn;
	delete [] ar;

	// delete single arrays
	delete [] Na;
	delete [] ac;
	delete [] cm;
	delete [] eulang1;
	delete [] eulang2;
	delete [] eulang3;

	// delete quaternion arrays
	delete [] q;
	delete [] qhalf;
	delete [] qnew;
	delete [] qdot;
}

/*
==================================

		  INITIALIZATION

==================================
*/

void rigidbody::get_file_header(string &rbstr) {
	int test,Ntest,i,j,l;
	double val, Ltmp, nasum;
	ifstream obj(rbstr.c_str());
	if (!obj.is_open()){
		cout << "Input file not opened, throwing error..." << endl;
		throw;
	}

	// determine what type of file (input or config)
	obj >> test;
	if (test == N){
		obj >> L[0] >> L[1] >> L[2];
		obj >> phi;

		// get rid of header
		char ctr;
		obj >> ctr;
		string header(" ");
		string equal_str(" ");
		getline(obj, header);
		getline(obj, equal_str);
		cout << "header = " << header << endl << endl;
		cout << equal_str << endl;
	}
	else{
		// get number of particles
		i = test;
		obj >> N;
		obj >> Ltmp;
		obj >> phi;

		// get L (already initialized through packing class);
		for (i = 0; i < NDIM; i++)
			L[i] = Ltmp;

		// get rid of header
		char ctr;
		obj >> ctr;
		string header(" ");
		string equal_str(" ");
		getline(obj, header);
		cout << "header = " << header << endl << endl;
	}	

	// loop over lines in data, get Na[i]
	Na = new int[N];
	i = 0;
	l = 0;
	Na[i] = 0;
	while (!obj.eof()) {
		obj >> i;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		obj >> val;
		// obj >> val;
		obj >> Na[i];
		cout << "l = " << l << "; i = " << i << "; Na[" << i << "] = " << Na[i] << endl;

		// increment l
		l++;

		// test to stop loop
		if (i == N-1){
			nasum = 0;
			for (j=0; j<N; j++)
				nasum += Na[j];
			if (l == nasum){
				cout << "all atoms found! breaking..." << endl;
				break;
			}
		}
	}

	// close file
	obj.close();
}

void rigidbody::initialize_particles() {
	// local variables
	int i, j;

	// initialize positions, atom info
	xW = new double**[N];
	xM = new double**[N];

	Inn = new double*[N];
	ar = new double*[N];

	for (i = 0; i < N; i++) {
		xW[i] = new double*[Na[i]];
		xM[i] = new double*[Na[i]];

		Inn[i] = new double[NDIM];
		ar[i] = new double[Na[i]];

		for (j = 0; j < Na[i]; j++) {
			xW[i][j] = new double[NDIM];
			xM[i][j] = new double[NDIM];
		}
	}

	// initialize quaternions
	q = new Quaternion[N];
	qhalf = new Quaternion[N];
	qnew = new Quaternion[N];
	qdot = new Quaternion[N];

	// initialize euler angles
	eulang1 = new double[N];
	eulang2 = new double[N];
	eulang3 = new double[N];

	// initialize other atomic info
	ac = new int[N];

	// initialize contact multiplicity matrix
	cm = new int[NC];
}

void rigidbody::initialize_dynamics() {
	// local variables
	int i, d;

	// allocate memory
	wW = new double*[N];
	wM = new double*[N];

	TqW = new double*[N];
	TqM = new double*[N];

	LW = new double*[N];
	LWhalf = new double*[N];
	LM = new double*[N];
	LMhalf = new double*[N];
	LMdot = new double*[N];

	// allocate second dimension
	for (i = 0; i < N; i++) {
		wW[i] = new double[NDIM];
		wM[i] = new double[NDIM];

		TqW[i] = new double[NDIM];
		TqM[i] = new double[NDIM];

		LW[i] = new double[NDIM];
		LWhalf[i] = new double[NDIM];
		LM[i] = new double[NDIM];
		LMhalf[i] = new double[NDIM];
		LMdot[i] = new double[NDIM];

		for (d = 0; d < NDIM; d++) {
			wW[i][d] = 0;
			wM[i][d] = 0;

			TqW[i][d] = 0;
			TqM[i][d] = 0;

			LW[i][d] = 0;
			LWhalf[i][d] = 0;
			LM[i][d] = 0;
			LMhalf[i][d] = 0;
			LMdot[i][d] = 0;
		}
	}
}

void rigidbody::read_in_info(string &rbstr) {
	int i,j,d,test;
	double rscale = 1;
	double val;
	ifstream obj(rbstr.c_str());

	// determine what type of file (input or config)
	obj >> test;
	if (test == N){
		obj >> val >> val >> val;
		obj >> val;

		// get rid of header
		char ctr;
		obj >> ctr;
		string header(" ");
		string equal_str(" ");
		getline(obj, header);
		getline(obj, equal_str);
		cout << "header = " << header << endl << endl;
		cout << equal_str << endl;

		// do NOT half radius
		rscale = 1;
	}
	else{
		// get number of particles
		obj >> val;
		obj >> val;
		obj >> val;

		// get rid of header
		char ctr;
		obj >> ctr;
		string header(" ");
		string equal_str(" ");
		getline(obj, header);
		cout << "header = " << header << endl << endl;

		// half radius
		rscale = 0.5;
	}

	// max residue radius
	double maxrad = 0;

	// loop over lines in data, get data
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			obj >> val;
			obj >> ar[i][j];
			ar[i][j] = rscale * ar[i][j];
			obj >> xW[i][j][0];
			obj >> xW[i][j][1];
			obj >> xW[i][j][2];
			obj >> r[i];
			r[i] = rscale * r[i];
			obj >> x[i][0];
			obj >> x[i][1];
			obj >> x[i][2];
			obj >> Inn[i][0];
			obj >> Inn[i][1];
			obj >> Inn[i][2];
			obj >> eulang1[i];
			obj >> eulang2[i];
			obj >> eulang3[i];
			obj >> m[i];
			obj >> val;

			if (r[i] > maxrad)
				maxrad = r[i];
		}
	}

	// make world frame position relative
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			for (d = 0; d < NDIM; d++)
				xW[i][j][d] = xW[i][j][d] - x[i][d];
		}
	}

	// close file
	obj.close();

	// if large enough box, update cell grid length
	if (NCL > 0) {
		init_rcut();
		update_cell_g();
	}
}

void rigidbody::initialize_quaternions() {
	// local variables
	int i, j;
	double qs, qx, qy, qz;

	// loop over particles
	for (i = 0; i < N; i++) {
		qs = cos(0.5 * eulang2[i]) * cos(0.5 * (eulang1[i] + eulang3[i]));
		qx = sin(0.5 * eulang2[i]) * cos(0.5 * (eulang1[i] - eulang3[i]));
		qy = sin(0.5 * eulang2[i]) * sin(0.5 * (eulang1[i] - eulang3[i]));
		qz = cos(0.5 * eulang2[i]) * sin(0.5 * (eulang1[i] + eulang3[i]));

		q[i].set_s(qs);
		q[i].set_x(qx);
		q[i].set_y(qy);
		q[i].set_z(qz);

		// normalize each quaternion
		q[i].normalize();
	}

	// update xM based on quaternions and xW
	pos_brot();
}

void rigidbody::rand_vel_init(double T1) {
	int i,d;
	double ek,vs;
	vector<double> pmean(NDIM,0.0);
	srand48(seed);

	// get initial velocities
	for (i=0; i<N; i++){
    	for (d=0; d<NDIM; d++){
        	v[i][d] = drand48();
        	pmean[d] += m[i]*v[i][d];        	
    	}
    }

    // get CoM momentum
    for (d=0; d<NDIM; d++)
		pmean[d] /= N;

	// subtract off CoM momentum
	for (i=0; i<N; i++){
    	for (d=0; d<NDIM; d++)
        	v[i][d] -= pmean[d]/m[i];
    }

    // rescale velocities
    rescale_velocities(T1);
}

// scramble particle positions and orientations
void rigidbody::scramble(){

	// local variables
	int i,j,d;

	// seed random number generator
	srand48(seed);

	// scramble particles
	for (i=0; i<N; i++){
		// scramble particle centers of mass
		for (d=0; d<NDIM; d++)
			x[i][d] = L[d]*drand48();

		// scramble orientations
		q[i].set_s(drand48());
		q[i].set_x(drand48());
		q[i].set_y(drand48());
		q[i].set_z(drand48());

		// normalize each quaternion
		q[i].normalize();
	}

	// update xM based on quaternions and xW
	pos_frot();
}


/*
==================================

		 SETTERS & GETTERS

==================================
*/

void rigidbody::reset_cm() {
	int i;

	for (i = 0; i < NC; i++)
		cm[i] = 0;
}

void rigidbody::update_phi() {
	int i, d;
	double msum, Lprod;

	msum = 0;
	Lprod = 1;
	for (i = 0; i < N; i++)
		msum += m[i];
	for (d = 0; d < NDIM; d++)
		Lprod *= L[d];

	phi = msum / Lprod;
}

void rigidbody::update_euler() {
	int i;
	double qs, qx, qy, qz;
	double se1,ce1,se2,ce2,se3,ce3;

	for (i = 0; i < N; i++) {
		// store quaternion values
		qs = q[i].get_s();
		qx = q[i].get_x();
		qy = q[i].get_y();
		qz = q[i].get_z();

		// get euler angle cosines and sines

		// polar angle
	    se2 = 2*sqrt((qx*qx + qy*qy)*(1-qx*qx-qy*qy));
	    ce2 = 1 - 2*(qx*qx+qy*qy);
	    
	    // azimuthal angle
	    se1 = 2*(qx*qz+qs*qy)/se2;
	    ce1 = 2*(qs*qx-qy*qz)/se2;
	    
	    // spin angle
	    se3 = 2*(qx*qz - qs*qy)/se2;
	    ce3 = 2*(qs*qx + qy*qz)/se2;

	    // update euler angles
	    eulang1[i] = atan2(se1,ce1);
	    eulang2[i] = atan2(se2,ce2);
	    eulang3[i] = atan2(se3,ce3);
	}
}

void rigidbody::update_quaternions() {
	// local variables
	int i, j;
	double qs, qx, qy, qz;

	// loop over particles
	for (i = 0; i < N; i++) {
		qs = cos(0.5 * eulang2[i]) * cos(0.5 * (eulang1[i] + eulang3[i]));
		qx = sin(0.5 * eulang2[i]) * cos(0.5 * (eulang1[i] - eulang3[i]));
		qy = sin(0.5 * eulang2[i]) * sin(0.5 * (eulang1[i] - eulang3[i]));
		qz = cos(0.5 * eulang2[i]) * sin(0.5 * (eulang1[i] + eulang3[i]));

		q[i].set_s(qs);
		q[i].set_x(qx);
		q[i].set_y(qy);
		q[i].set_z(qz);

		// normalize each quaternion
		q[i].normalize();
	}

	// update xW based on quaternions and xM and quaternions
	pos_frot();
}

int rigidbody::get_ac_sum() {
	int i;
	int val = 0;
	for (i = 0; i < N; i++)
		val += ac[i];
	return val;
}

double rigidbody::get_atomic_distance(int i, int j, int ai, int aj, double aij[]) {
	int d;
	double dr, h, ap1, ap2;

	h = 0;
	dr = 0;
	for (d = 0; d < NDIM; d++) {
		// get position in W coordinate, not relative
		ap1 = xW[i][ai][d] + x[i][d];
		ap2 = xW[j][aj][d] + x[j][d];

		// get distance
		dr = ap2 - ap1;
		dr -= L[d] * round(dr / L[d]);

		// save distance for force calc
		aij[d] = dr;
		h += dr * dr;
	}

	h = sqrt(h);
	return h;
}

double rigidbody::get_Natot() {
	int i;
	double sum = 0;

	for (i = 0; i < N; i++)
		sum += Na[i];

	return sum;
}

double rigidbody::get_LWX() {
	int i;
	double lwsum = 0;

	for (i = 0; i < N; i++)
		lwsum += LW[i][0];

	return lwsum;
}

double rigidbody::get_LWY() {
	int i;
	double lwsum = 0;

	for (i = 0; i < N; i++)
		lwsum += LW[i][1];

	return lwsum;
}

double rigidbody::get_LWZ() {
	int i;
	double lwsum = 0;

	for (i = 0; i < N; i++)
		lwsum += LW[i][2];

	return lwsum;
}


/*

==================================

			MD/PACKING
			PROTOCOLS

==================================

*/

void rigidbody::free_fire(double tmp0, double Utol, double tend, int nnu) {
	int t;

	// initialize velocities
	rand_vel_init(tmp0);

	// setup dtmax for FIRE
	set_dtmax(10.0);

	// get number of time steps
	int NT;
	NT = round(tend / dt);
	dt = tend / NT;

	// update nearest neighbor update
	nnupdate = nnu;

	// if energy output open, output dt
	enobj << dt << endl;

	cout << "===== BEGINNING TIME LOOP ======" << endl;
	cout << "NT = " << NT << endl;
	cout << "dt = " << dt << endl;
	cout << "tend = " << tend << endl;

	for (t = 0; t < NT; t++) {
		// if NLCL, update
		if (NCL > 0 && t % nnupdate == 0) {
			update_cell();
			update_neighborlist();
		}

		// advance quaternions, positions
		verlet_first();

		// update forces
		force_update();

		if (U < Utol)
			break;

		// include fire relaxation
		rb_fire();

		// advance angular momentum
		verlet_second();

		// output information
		if (t % plotskip == 0) {
			monitor_header(t);
			rigidbody_md_monitor();
		}
	}
}

void rigidbody::rb_jamming_finder(double tmp0, int NT, double dphi, double Utol, double Ktol) {
	// local variables
	int t, kr, check_rattlers;
	double dphi0, phiL, phiH;

	// initialize variables
	kr = 0;
	dphi0 = dphi;
	phiL = -1;
	phiH = -1;
	check_rattlers = 0;

	// constant energy checking
	double Uold, dU, dUtol;
	int epc, epcN, epconst;
	epconst = 0;
	Uold = 0;
	dU = 0;
	dUtol = 1e-8;
	epc = 0;
	epcN = 5e2;

	// make nnupdate larger
	nnupdate = 20;

	// initialize velocities
	rand_vel_init(tmp0);

	// setup dtmax for FIRE
	set_dtmax(10.0);

	cout << "===== BEGINNING TIME LOOP ======" << endl;
	cout << "NT = " << NT << endl;
	cout << "dt = " << dt << endl;
	cout << "================================" << endl << endl;

	for (t = 0; t < NT; t++) {
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
		if (check_rattlers) {
			kr = 0;
			nr = rmv_rattlers(kr);
		}
		else
			nr = 0;

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

		if (epconst == 1)
			check_rattlers = 1;

		// run root search routine
		rb_root_search(phiH, phiL, check_rattlers, epconst, nr, dphi0, Ktol, Utol, t);

		// output information
		if (t % plotskip == 0) {
			monitor_header(t);
			rigidbody_md_monitor();
			cout << "** ROOT SEARCH: " << endl;
			cout << "phi = " << phi << endl;
			cout << "phiL = " << phiL << endl;
			cout << "phiH = " << phiH << endl;
			cout << "epconst = " << epconst << endl;
			cout << "check_rattlers = " << check_rattlers << endl;
			cout << "nr = " << nr << endl;
			cout << endl;
			cout << endl;
		}

		// break if jamming found
		if (isjammed == 1)
			break;
		else if (isjammed == -1){
			cout << "isjammed was found to be -1, ending..." << endl;
			break;
		}
	}
}

void rigidbody::rb_anneal(double tmp0, int NT, int fskip, double phimin, double dphi, double Utol, double Ktol){
	// local variables
	int t, kr, check_rattlers, kskip, anneal, checkphi;
	double dphi0, phiL, phiH;

	// initialize variables
	kr = 0;
	dphi0 = dphi;
	phiL = -1;
	phiH = -1;
	check_rattlers = 0;

	// constant energy checking
	double Uold, dU, dUtol;
	int epc, epcN, epconst;
	epconst = 0;
	Uold = 0;
	dU = 0;
	dUtol = 1e-12;
	epc = 0;
	epcN = 5e3;

	// annealing procedure: skip fire/growth after fskip annealing steps
	kskip = 0;
	anneal = 0;
	checkphi = 1;

	// initialize velocities
	rand_vel_init(tmp0);

	// setup dtmax for FIRE
	set_dtmax(10.0);

	cout << "===== BEGINNING TIME LOOP ======" << endl;
	cout << "NT = " << NT << endl;
	cout << "dt = " << dt << endl;
	cout << "================================" << endl << endl;

	for (t = 0; t < NT; t++) {
		// update nearest neighbor lists if applicable
		if (t % nnupdate == 0 && NCL > -1){
			cout << "^ ";
			update_nlcl(t);
		}

		// only begin annealing after phi > phimin
		if (phi > phimin && checkphi == 1){
			anneal = 1;
			checkphi = 0;
		}

		// advance quaternions, positions
		verlet_first();

		// update forces
		force_update();

		// include fire relaxation (if not annealing)
		if (anneal == 0)
			rb_fire();
		else
			rescale_velocities(tmp0);
			

		// advance angular momentum
		verlet_second();		

		// check for rattlers
		if (check_rattlers) {
			kr = 0;
			nr = rmv_rattlers(kr);
		}
		else
			nr = 0;

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

		if (epconst == 1)
			check_rattlers = 1;

		// run root search routine (if not annealing)
		if (anneal == 0){
			rb_root_search(phiH, phiL, check_rattlers, epconst, nr, dphi0, Ktol, Utol, t);

			// if still growing, and energy relaxed, go back to annealing
			if (phiH < 0 && U < Utol && checkphi ==  0){
				anneal = 1;

				// give random velocity kick
				rand_vel_init(tmp0);
			}
		}
		else{
			kskip++;
			if (kskip == fskip){
				anneal = 0;
				kskip = 0;
			}
		}

		// output information
		if (t % plotskip == 0) {
			monitor_header(t);
			rigidbody_md_monitor();
			cout << "** ROOT SEARCH: " << endl;
			cout << "phi = " << phi << endl;
			cout << "phiL = " << phiL << endl;
			cout << "phiH = " << phiH << endl;
			cout << "epconst = " << epconst << endl;
			cout << "check_rattlers = " << check_rattlers << endl;
			cout << "nr = " << nr << endl;
			cout << endl;
			cout << endl;
		}

		// break if jamming found
		if (isjammed == 1)
			break;
		else if (isjammed == -1){
			cout << "isjammed was found to be -1, ending..." << endl;
			break;
		}
	}
}












void rigidbody::rb_jamming_precise(double tphiold, double tmp0, int NT, double Utol, double Ktol) {
	// local variables
	int s, kr, check_rattlers, nr;
	double dphi0, phiold, phinew, tphinew1, tphinew2, tphinew, Uold, Unew;
	double y = pow(10,-1.0/10.0);

	// initialize velocities
	rand_vel_init(tmp0);

	// initialize variables
	nr = 0;
	phiold = 2*tphiold;

	// get initial U
	get_U(Ktol,nr);
	Uold = U;

	// update phinew
	phinew = tphiold + (phiold-tphiold)*y;	

	// update system
	rb_scale(phinew);

	// get final U
	get_U(Ktol,nr);
	Unew = U;

	// make nnupdate larger
	nnupdate = 20;

	// get new guess for phiJ
	tphinew1 = phinew - phiold*sqrt(Unew/Uold);
	tphinew2 = 1 - sqrt(Unew/Uold);
	tphinew = tphinew1/tphinew2;

	// update phis for next run
	phiold = phinew;
	tphiold = tphinew;
	Uold = Unew;

	// loop over time
	isjammed = 0;
	for (s=0; s<NT; s++){
		// update phinew
		phinew = tphiold + (phiold-tphiold)*y;			

		// update system
		rb_scale(phinew);

		// get updated U
		get_U(Ktol,nr);
		Unew = U;

		// output to energy file if open
		if (enobj.is_open()) {
			cout << "Printing ENERGY" << endl;
			enobj << setw(12) << U;
			enobj << setw(12) << K;
			enobj << setw(12) << Krot;
			enobj << setw(12) << U + K;
			enobj << endl;
		}

		// get new guess for phiJ
		tphinew1 = phinew - phiold*sqrt(Unew/Uold);
		tphinew2 = 1 - sqrt(Unew/Uold);
		tphinew = tphinew1/tphinew2;

		cout << endl << endl;
		cout << "###################" << endl;
		cout << "## s  = " << s << endl;
		cout << "## phinew = " << setw(12) << phinew << "; phiold = " << setw(12) << phiold << endl;
		cout << "## tphinew = " << setw(12) << tphinew << "; tphiold = " << setw(12) << tphiold << endl;
		cout << "## new U = " << setw(12) << Unew << "; Uold = " << setw(12) << Uold << endl;
		cout << "## nr = " << nr << endl;
		cout << "## niso = " << DOF*(N-nr) - NDIM + 1 << endl;
		cout << "## acsum = " << 0.5*get_ac_sum() << endl;
		cout << "## pcsum = " << get_c_sum() << endl;				
		cout << "###################" << endl;
		cout << endl << endl;


		// check for convergence
		if (U < Utol){
			cout << "Jammed packing found!" << endl;
			cout << endl;
			cout << endl;
			cout << "Found jammed state!" << endl;
			cout << "Writing config after s = " << s << " steps..." << endl;
			cout << "Final U = " << U << endl;
			cout << "Final K = " << K << endl;
			cout << "Final phi = " << setprecision(6) << phi << endl;
			cout << "Final acsum = " << 0.5*get_ac_sum() << endl;
			cout << "Final pcsum = " << get_c_sum() << endl;
			cout << "Final niso max = " << DOF*N - NDIM + 1 << endl;
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

			cout << endl;
			cout << endl;
			isjammed = 1;
			break;
		}		

		// update phis for next run
		phiold = phinew;
		tphiold = tphinew;
		Uold = Unew;
	}

	if (s == NT)
		cout << "Jammed state was not found in NT..." << endl;
}

void rigidbody::rb_jamming_easy(double tmp0, int NT, double dphi, double Utol, double Ktol) {
	// local variables
	int t, kr, check_rattlers;
	double dphi0, phiL, phiH;
	bool min = 0;

	// initialize variables
	kr = 0;
	dphi0 = dphi;
	phiL = -1;
	phiH = -1;
	check_rattlers = 0;

	// constant energy checking
	double Uold, dU, dUtol;
	int epc, epcN, epconst;
	epconst = 0;
	Uold = 0;
	dU = 0;
	dUtol = 1e-8;
	epc = 0;
	epcN = 5e3;

	// make nnupdate larger
	nnupdate = 20;

	// initialize velocities
	rand_vel_init(tmp0);

	// setup dtmax for FIRE
	set_dtmax(10.0);

	cout << "===== BEGINNING TIME LOOP ======" << endl;
	cout << "NT = " << NT << endl;
	cout << "dt = " << dt << endl;
	cout << "================================" << endl << endl;

	for (t = 0; t < NT; t++) {
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
		if (check_rattlers) {
			kr = 0;
			nr = rmv_rattlers(kr);
		}
		else
			nr = 0;

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

		if (epconst == 1)
			check_rattlers = 1;

		// run root search routine
		rb_easy(phiH, phiL, check_rattlers, epconst, nr, dphi0, Ktol, Utol, t, min);

		// output information
		if (t % plotskip == 0) {
			monitor_header(t);
			rigidbody_md_monitor();
			cout << "** ROOT SEARCH: " << endl;
			cout << "phi = " << phi << endl;
			cout << "phiL = " << phiL << endl;
			cout << "phiH = " << phiH << endl;
			cout << "epconst = " << epconst << endl;
			cout << "check_rattlers = " << check_rattlers << endl;
			cout << "nr = " << nr << endl;
			cout << endl;
			cout << endl;
		}

		// break if jamming found
		if (isjammed == 1)
			break;
	}
}

void rigidbody::get_U(double Ktol, int& nr){
	int t,NT,kr;	
	NT = 5e6;


	// constant energy checking
	double Uold, dU, dUtol;
	int epc, epcN, epconst;
	epconst = 0;
	Uold = 0;
	dU = 0;
	dUtol = 1e-8;
	epc = 0;
	epcN = 1e2;
	
	// Minimize energy using FIRE
	cout << "** GETTING Unew FOR phi = " << phi << endl;
	for (t = 0; t < NT; t++) {
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
			monitor_header(t);
			rigidbody_md_monitor();
			cout << "U = " << U << endl;
			cout << "Uold = " << Uold << endl;
			cout << "dU = " << dU << endl;
			cout << "epconst = " << epconst << endl;
			cout << "nr = " << nr << endl;
			cout << endl;
			cout << endl;
		}

		if (epconst == 1 && K < Ktol)
			break;
	}

	if (t==NT)
		cout << "** COULD NOT FIND ENERGY MINIMUM IN NT', ENDING AND SETTING UNEW = U..." << endl;
}


/*

==================================

			MD/PACKING
		      CODE

==================================

*/


void rigidbody::verlet_first() {
	// update translational pos first, to get kinetic energy rolling
	pos_update();

	// update quaternions second, and include Krot update in euler_q
	q_step();

	// reset rotational variables
	int i, d;
	for (i = 0; i < N; i++) {
		// reset atomic contacts
		ac[i] = 0;
		for (d = 0; d < NDIM; d++) {
			// reset torques
			TqW[i][d] = 0;
		}
	}
}

void rigidbody::verlet_second(bool neglect_rattlers) {
	int i, d;
	double vel, anew;

	for (i = 0; i < N; i++) {

		// if input exists and is true, then do not update rattler velocities
		if (neglect_rattlers){
			if (pc[i] == 0)
				continue;
		}

		for (d = 0; d < NDIM; d++) {
			// update angular momentum
			LW[i][d] = LWhalf[i][d] + 0.5 * dt * TqW[i][d];

			// update velocties
			vel = v[i][d];
			anew = F[i][d] / m[i];
			vel += 0.5 * dt * (anew + aold[i][d]);
			v[i][d] = vel;
			aold[i][d] = anew;
		}
	}
}

// quaternion-MD Step
void rigidbody::q_step() {
	// 1. Rotate W to M frame
	rotation_W2M();

	// 2. Euler step for quaternions
	euler_q();

	// 3. Update atomic positions in World frame
	pos_frot();
}

void rigidbody::rotation_W2M() {
	// local variables
	int i, d;

	// angular momentum and torque objects
	Quaternion lwtmp;
	Quaternion twtmp;

	// placeholder objects
	Quaternion q1;
	Quaternion q2;

	for (i = 0; i < N; i++) {
		// setup temporary quaternions
		lwtmp.set_s(0.0);
		lwtmp.set_x(LW[i][0]);
		lwtmp.set_y(LW[i][1]);
		lwtmp.set_z(LW[i][2]);

		twtmp.set_s(0.0);
		twtmp.set_x(TqW[i][0]);
		twtmp.set_y(TqW[i][1]);
		twtmp.set_z(TqW[i][2]);

		// rotate with quaternion algebra
		q1 = lwtmp % q[i];
		q2 *= q[i];
		lwtmp = q2 % q1;

		q1 = twtmp % q[i];
		q2 *= q[i];
		twtmp = q2 % q1;

		// update molecule frame
		LM[i][0] = lwtmp.get_x();
		LM[i][1] = lwtmp.get_y();
		LM[i][2] = lwtmp.get_z();

		TqM[i][0] = twtmp.get_x();
		TqM[i][1] = twtmp.get_y();
		TqM[i][2] = twtmp.get_z();
	}
}

void rigidbody::euler_q() {
	// local variables
	int i, d;
	double cx, cy, cz;
	Quaternion qtmp;
	Quaternion q1;
	Quaternion q2;
	Quaternion qold;

	// tolerance
	double qep = 1e-8;
	double qcheck = 10 * qep;
	int k = 0;
	int kmax = 1e4;

	// restart Krot
	Krot = 0;

	// update w, Ldot, advance to half step
	for (i = 0; i < N; i++) {
		// angular velocity
		for (d = 0; d < NDIM; d++)
			wM[i][d] = LM[i][d] / Inn[i][d];

		// L time derivative
		cx = wM[i][1] * LM[i][2] - wM[i][2] * LM[i][1];
		cy = -wM[i][0] * LM[i][2] + wM[i][2] * LM[i][0];
		cz = wM[i][0] * LM[i][1] - wM[i][1] * LM[i][0];

		LMdot[i][0] = TqM[i][0] - cx;
		LMdot[i][1] = TqM[i][1] - cy;
		LMdot[i][2] = TqM[i][2] - cz;

		// half euler step for LM
		for (d = 0; d < NDIM; d++)
			LMhalf[i][d] = LM[i][d] + 0.5 * dt * LMdot[i][d];

		// approx qdot
		qtmp.set_s(0.0);
		qtmp.set_x(0.5 * (LMhalf[i][0] / Inn[i][0]));
		qtmp.set_y(0.5 * (LMhalf[i][1] / Inn[i][1]));
		qtmp.set_z(0.5 * (LMhalf[i][2] / Inn[i][2]));
		qdot[i] = q[i] % qtmp;

		// half euler step for q
		qhalf[i] = q[i] + qdot[i] * (0.5 * dt);
		qhalf[i].normalize();

		// half euler step for worlf frame L (LW)
		for (d = 0; d < NDIM; d++)
			LWhalf[i][d] = LW[i][d] + 0.5 * dt * TqW[i][d];

		// determine qdot self-consistently
		while (qcheck > qep && k < kmax) {
			// update k
			k++;

			// store old qhalf
			qold = qhalf[i];

			// rotate LWhalf to LMhalf
			qtmp.set_s(0.0);
			qtmp.set_x(LWhalf[i][0]);
			qtmp.set_y(LWhalf[i][1]);
			qtmp.set_z(LWhalf[i][2]);

			q1 = qtmp % qhalf[i];
			q2 *= qhalf[i];
			qtmp = q2 % q1;

			LMhalf[i][0] = qtmp.get_x();
			LMhalf[i][1] = qtmp.get_y();
			LMhalf[i][2] = qtmp.get_z();

			// update angular velocity
			for (d = 0; d < NDIM; d++)
				wM[i][d] = LMhalf[i][d] / Inn[i][d];

			// approx qdot
			qtmp.set_s(0.0);
			qtmp.set_x(0.5 * wM[i][0]);
			qtmp.set_y(0.5 * wM[i][1]);
			qtmp.set_z(0.5 * wM[i][2]);
			qdot[i] = qhalf[i] % qtmp;

			// half euler step for q
			qhalf[i] = q[i] + qdot[i] * (0.5 * dt);

			// update qcheck
			qold = qhalf[i] - qold;
			qcheck = qold.get_norm();
		}

		// take full euler step for q
		q[i] = q[i] + qdot[i] * dt;

		// enforce quaternion normalization
		q[i].normalize();

		// update rotational kinetic energy
		for (d = 0; d < NDIM; d++)
			Krot += 0.5 * Inn[i][d] * pow(wM[i][d], 2);
	}
	K += Krot;
}

void rigidbody::pos_frot() {
	int i, j, d;
	Quaternion q1;
	Quaternion q2;
	Quaternion xtmp;

	// rotate every particle position in M frame
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			// rotate rel M frame to rel W frame
			xtmp.set_s(0.0);
			xtmp.set_x(xM[i][j][0]);
			xtmp.set_y(xM[i][j][1]);
			xtmp.set_z(xM[i][j][2]);

			q1 *= q[i];
			q2 = xtmp % q1;
			xtmp = q[i] % q2;

			xW[i][j][0] = xtmp.get_x();
			xW[i][j][1] = xtmp.get_y();
			xW[i][j][2] = xtmp.get_z();
		}
	}
}

void rigidbody::pos_brot() {
	int i, j, d;
	Quaternion q1;
	Quaternion q2;
	Quaternion xtmp;

	// rotate every particle position in w frame
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			// rotate rel W frame to rel M frame
			xtmp.set_s(0.0);
			xtmp.set_x(xW[i][j][0]);
			xtmp.set_y(xW[i][j][1]);
			xtmp.set_z(xW[i][j][2]);

			// rotate xW to xM through backward rotation
			q1 = xtmp % q[i];
			q2 *= q[i];
			xtmp = q2 % q1;

			xM[i][j][0] = xtmp.get_x();
			xM[i][j][1] = xtmp.get_y();
			xM[i][j][2] = xtmp.get_z();
		}
	}
}

void rigidbody::rotate_single_particle(int i, int angle, double dtheta){

	// first make sure euler angles are up to date
	update_euler();

	// rotate given angle by dtheta
	if (angle == 1)
		eulang1[i] += dtheta;
	else if (angle == 2)
		eulang2[i] += dtheta;
	else if (angle == 3)
		eulang3[i] += dtheta;
	else{
		cout << "not applicable angle = " << angle << "..." << endl;
		exit(1);
	}

	// update quaternions accordingly
	update_quaternions();
}

void rigidbody::rotate_single_particle_xyz(int i, int axis, double dtheta){
	// local variables
	double ci,si; 	// cosines of angle
	double vj, vk; 	// rotated angles
	int j,k,a;		// other directions

	double vnorm1, vnorm2;

	// cosines and sines
	ci = cos(dtheta);
	si = sin(dtheta);
	cout << "dtheta = " << dtheta << endl;
	cout << "ci = " << ci << endl;
	cout << "si = " << si << endl;

	// test axis
	if (axis > 2 || axis < 0){
		cout << "	** ERROR: incorrect axis passed as argument (" << axis << "), which is out of bounds. Ending." << endl;
		exit(1);
	}

	// perform rotation
	j = (axis + 1) % NDIM;
	k = (axis + 2) % NDIM;

	// get updated vectors
	for (a=0; a<Na[i]; a++){
		// get new components
		vj = ci*xW[i][a][j] - si*xW[i][a][k];
		vk = si*xW[i][a][j] + ci*xW[i][a][k];

		// update positions
		xW[i][a][j] = vj;
		xW[i][a][k] = vk;
	}

	
}

void rigidbody::force_update() {
	int i, j, jj, ai, aj, d, cind, M, pc_found;
	double sij, rij, dx, da;
	double qix, qiy, qiz, qjx, qjy, qjz;
	double tix, tiy, tiz, tjx, tjy, tjz;
	double fix, fiy, fiz;
	double Rij[NDIM];
	double aij[NDIM];
	double fij[NDIM];

	// reset U, LCON
	U = 0;
	LCON = 0;
	reset_c();
	reset_cm();

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
				dx = get_distance(i, j);

				if (dx < sij) {
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
									fij[d] = (hs(rij, da)) * aij[d];
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
								U += (ep / 2) * pow(1 - da / rij, 2);
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

					// check local angular momentum conservation (LCON should be zero in all 3 directions)
					LCON += tix + tjx + tiy + tjy + tiz + tjz;
					LCON += -Rij[1] * fiz + Rij[2] * fiy;
					LCON += Rij[0] * fiz - Rij[2] * fix;
					LCON += -Rij[0] * fiy + Rij[1] * fix;
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
				dx = get_distance(i, j);

				// if true, residues close by, check atomic overlaps
				if (dx < sij) {
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
									fij[d] = (hs(rij, da)) * aij[d];
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
								U += (ep / 2) * pow(1 - da / rij, 2);
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

					// check local angular momentum conservation (LCON should be zero in all 3 directions)
					LCON += tix + tjx + tiy + tjy + tiz + tjz;
					LCON += -Rij[1] * fiz + Rij[2] * fiy;
					LCON += Rij[0] * fiz - Rij[2] * fix;
					LCON += -Rij[0] * fiy + Rij[1] * fix;
				}
			}
		}
	}
}

// remove particles with unconstranted dof
int rigidbody::rmv_rattlers(int& krcrs) {
	int i, j, ci, cj, r, nr, nm;

	// monitor recursion depth
	krcrs++;

	// number of rattlers
	nr = 0;

	// number of "marginal" rattlers to be removed
	nm = 0;

	// loop over rows, eliminate contacts to rattlers
	for (i = 0; i < N; i++) {
		// get number of contacts
		r = ac[i];

		// remove from network if r <= DOF, delete contacts
		if (r <= DOF) {
			// increment # of rattlers
			nr++;

			// alter contact vectors
			ac[i] = 0;
			pc[i] = 0;

			// if in contact, remove contacts
			if (r > 0) {
				nm++;
				for (j = 0; j < N; j++) {
					if (j < i) {
						ci = N * j + i - ((j + 1) * (j + 2)) / 2;
						if (c[ci] > 0) {
							pc[j]--;
							ac[j] -= cm[ci];
							c[ci] = 0;
							cm[ci] = 0;
						}
					}
					else if (j > i) {
						cj = N * i + j - ((i + 1) * (i + 2)) / 2;	// mapping from matrix space to sub matrix space
						if (c[cj] > 0) {
							pc[j]--;
							ac[j] -= cm[cj];
							c[cj] = 0;
							cm[cj] = 0;
						}
					}
					else
						continue;
				}
			}
		}
	}

	if (krcrs > 100) {
		cout << "max recursive depth reached, be wary of output" << endl;
		return -1;
	}

	// recursively check once rattlers are removed
	if (nm == 0)
		return nr;
	else
		return rmv_rattlers(krcrs);
}

int rigidbody::id_rattlers() {
	int i, dof_tol, nr;
	dof_tol = DOF;

	nr = 0;
	for (i = 0; i < N; i++) {
		if (ac[i] < dof_tol)
			nr++;
	}

	return nr;
}





/*

==================================

			GROWTH/EM
		      CODE

==================================

*/

// GROWTH methods
void rigidbody::rb_scale(double phinew) {
	int i, j, d;
	double s, invdim;

	// get scale parameter
	invdim = (double)1 / NDIM;
	s = pow(phinew / phi, invdim);

	// loop over system, scale everything
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			// scale positions (need to subtract off position increase)
			for (d = 0; d < NDIM; d++) {
				xM[i][j][d] *= s;
				xW[i][j][d] *= s;
			}

			// scale radii
			ar[i][j] *= s;
		}
		// scale residue mass
		m[i] *= pow(s, NDIM);

		// scale moment of inertia tensor
		for (d = 0; d < NDIM; d++)
			Inn[i][d] *= pow(s, NDIM + 2);

		// scale shell radii
		r[i] *= s;
	}
	ep *= pow(s, 2);
	dt *= pow(s, 0.5 * (NDIM-2));
	dtmax *= pow(s, 0.5 * (NDIM-2));
	update_phi();

	// if NLCL engaged, scale rcut
	if (NCL > -1) {
		scale_rcut(s);
		if (NCL > 3)
			update_cell();

		update_neighborlist();
	}
}

// rescale system by length len
void rigidbody::rb_rescale(double len) {
	int i,j,d;

	for (d=0; d<NDIM; d++){
		L[d] /= len;
	}

	// loop over system, scale everything
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			// scale positions (need to subtract off position increase)
			for (d = 0; d < NDIM; d++) {
				xM[i][j][d] /= len;
				xW[i][j][d] /= len;
			}

			// scale radii
			ar[i][j] /= len;
		}
		// scale residue mass
		m[i] /= pow(len, NDIM);

		// scale moment of inertia tensor
		for (d = 0; d < NDIM; d++)
			Inn[i][d] /= pow(len, NDIM+2);

		// scale shell radii
		r[i] /= len;

		// scale center of mass position
		for (d=0; d<NDIM; d++)
			x[i][d] /= len;
	}

	// if NLCL engaged, scale rcut
	if (NCL > -1) {
		scale_rcut(1.0/len);
		if (NCL > 3)
			update_cell();

		update_neighborlist();
	}
}

void rigidbody::rescale_velocities(double T1){
	// local variables
	int i,d;
	double T0,tmpScale;
	Quaternion qtmp, q1, q2;	

	// get current kinetic energy
	T0 = 0.0;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			T0 += 0.5*m[i]*v[i][d]*v[i][d] + 0.5*Inn[i][d]*wM[i][d]*wM[i][d];
	}

	// get scale for new kinetic energy
	tmpScale = sqrt(T1/T0);

	// rescale velocities and angular momenta
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++){
			v[i][d] *= tmpScale;
			wM[i][d] *= tmpScale;
			LM[i][d] = Inn[i][d]*wM[i][d];

			// rotate LM into LW
			qtmp.set_x(LM[i][0]);
			qtmp.set_y(LM[i][1]);
			qtmp.set_z(LM[i][2]);

			q1 *= q[i];
			q2 = qtmp % q1;
			qtmp = q[i] % q2;

			// update molecule frame
			LWhalf[i][0] = qtmp.get_x();
			LWhalf[i][1] = qtmp.get_y();
			LWhalf[i][2] = qtmp.get_z();
		}			
	}
}

// RIGIDBODY FIRE
void rigidbody::rb_fire() {
	int i, j, d;
	double P,Pi,Pimin;
	double vstarnrm = 0;
	double wstarnrm = 0;
	double fstarnrm = 0;
	double tstarnrm = 0;
	int Nmin = 50;
	Quaternion qtmp, q1, q2;	

	// calculate P
	P = 0;
	Pi = 0;
	Pimin = 1e8;
	for (i = 0; i < N; i++) {
		// rotate torques into M frame
		qtmp.set_x(TqW[i][0]);
		qtmp.set_y(TqW[i][1]);
		qtmp.set_z(TqW[i][2]);

		q1 = qtmp % q[i];
		q2 *= q[i];
		qtmp = q2 % q1;

		TqM[i][0] = qtmp.get_x();
		TqM[i][1] = qtmp.get_y();
		TqM[i][2] = qtmp.get_z();

		for (d = 0; d < NDIM; d++) {
			Pi = F[i][d] * v[i][d] + TqM[i][d] * wM[i][d];
			P += Pi;
			vstarnrm += v[i][d] * v[i][d];
			wstarnrm += wM[i][d] * wM[i][d];
			fstarnrm += F[i][d] * F[i][d];
			tstarnrm += TqM[i][d] * TqM[i][d];
		}
	}

	vstarnrm = sqrt(vstarnrm);
	wstarnrm = sqrt(wstarnrm);
	fstarnrm = sqrt(fstarnrm);
	tstarnrm = sqrt(tstarnrm);

	// update v if forces acting
	if (fstarnrm > 0 && tstarnrm > 0) {
		for (i = 0; i < N; i++) {
			for (d = 0; d < NDIM; d++) {
				v[i][d] = (1 - alpha) * v[i][d] + alpha * (F[i][d] / fstarnrm) * vstarnrm;
				wM[i][d] = (1 - alpha) * wM[i][d] + alpha * (TqM[i][d] / tstarnrm) * wstarnrm;
				LM[i][d] = Inn[i][d] * wM[i][d];
			}

			// rotate LM to LW
			qtmp.set_x(LM[i][0]);
			qtmp.set_y(LM[i][1]);
			qtmp.set_z(LM[i][2]);

			q1 *= q[i];
			q2 = qtmp % q1;
			qtmp = q[i] % q2;

			// update molecule frame
			LWhalf[i][0] = qtmp.get_x();
			LWhalf[i][1] = qtmp.get_y();
			LWhalf[i][2] = qtmp.get_z();
		}
	}

	// now update alphas for P
	if (P >= 0 && np > Nmin){

		// increase dt
		if (dt * finc < dtmax)
			dt *= finc;
		else
			dt = dtmax;

		// decrease alpha
		alpha *= falpha;

		np++;
	}
	else if (P < 0) {
		// reset K to measure based on new info
		K = 0;

		// decrease time step
		if (dt * fdec > 1e-2 * dtmax)
			dt *= fdec;
		else
			dt = 1e-2 * dtmax;

		// set global velocity vector to zero
		for (i = 0; i < N; i++) {
			for (d = 0; d < NDIM; d++) {
				K += 0.5 * m[i] * v[i][d] * v[i][d] + 0.5 * Inn[i][d] * wM[i][d] * wM[i][d];
				v[i][d] = 0;
				wM[i][d] = 0;
				LW[i][d] = 0;
				LWhalf[i][d] = 0;
				LM[i][d] = 0;
				LMhalf[i][d] = 0;
			}
			qdot[i].set_s(0);
			qdot[i].set_x(0);
			qdot[i].set_y(0);
			qdot[i].set_z(0);
		}

		// set alpha -> alphaStart
		alpha = alpha0;

		// set np -> 0
		np = 0;
	}
	else if (P >= 0 && np <= Nmin) 
		np++;
}

void rigidbody::rb_root_search(double& phiH, double& phiL, int& check_rattlers, int epconst, int nr, double dphi0, double Ktol, double &Utol, int t) {
	/*
		GROWTH ALGORITHM

		1. If U < Utol grow
		2. If U > 2*Utol & K < Ktol, stop, shrink by dphi0/2
			a. If U < Utol & K < Ktol, relaxed, so grow again
			b. If U > 2*Utol & K < Ktol, still overcompressed, so shrink by dphi0/2
		3. Jammed when K < Ktol & Utol < U < 2Utol
	*/

	int nbb, niso, pcsum, acsum;
	bool gr, oc, uc, jammed;
	double dphi = dphi0;

	gr = 0;
	oc = 0;
	uc = 0;
	jammed = 0;

	nbb = N - nr;
	niso = DOF * nbb - NDIM + 1;
	acsum = 0.5*get_ac_sum();
	pcsum = get_c_sum();

	gr = (U < Utol);
	oc = (U > 2 * Utol && K < Ktol && epconst == 1 && acsum > 0);
	uc = (U < Utol && epconst == 1);
	jammed = ( (U > Utol && U < 2 * Utol && K < Ktol && epconst == 1 && acsum > 0) );


	if (phiH < 0) {
		if (gr) {
			check_rattlers = 0;
			rb_scale(phi+dphi);
		}
		else if (oc && epconst == 1) {
			phiH = phi;
			dphi = -dphi0;
			check_rattlers = 1;

			cout << endl;
			cout << "phiH 1st set at nt = " << t << endl;
			rb_scale(phi+dphi);

			// if NLCL, change update check (particles don't move, don't need to check as often)
			if (NCL > -1)
				nnupdate *= 50;
		}
	}
	else {
		if (phiL < 0) {

			// if still overcompressed, decrease again
			if (oc && epconst == 1) {
				phiH = phi;
				dphi = -0.5*dphi0;				
				rb_scale(phi+dphi);
			}

			// if undercompressed, set phiL, root search
			if (uc) {
				phiL = phi;
				dphi = 0.5 * (phiH + phiL) - phi;
				rb_scale(phi+dphi);
			}

			if (jammed) {
				phiL = 0.99 * phi;
				dphi = 0.5 * (phiH + phiL) - phi;
				rb_scale(phi+dphi);
			}

		}
		else {

			// if overcompressed, root search down
			if (oc) {
				phiH = phi;
				dphi = 0.5 * (phiH + phiL) - phi;
				rb_scale(phi+dphi);
			}

			// if undercompressed, root search up
			if (uc) {
				phiL = phi;
				dphi = 0.5 * (phiH + phiL) - phi;
				rb_scale(phi+dphi);
			}

			// if jammed, end!
			if (jammed) {
				cout << endl;
				cout << endl;
				cout << "Found jammed state!" << endl;
				cout << "Writing config at t = " << t*dt << endl;
				cout << "Writing config at nt = " << t << endl;
				cout << "Final U = " << U << endl;
				cout << "Final K = " << K << endl;
				cout << "Final phi = " << setprecision(6) << phi << endl;
				cout << "Final pcsum = " << pcsum << endl;
				cout << "Final acsum = " << acsum << endl;
				cout << "Final niso = " << niso << endl;
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

				cout << endl;
				cout << endl;
				isjammed = 1;
			}
		}
	}

	// test for stalled growth
	if (abs(dphi) < 1e-14) {
		cout << " ** broken jamming found, no state found..." << endl;
		cout << " ** setting isjammed to -1, ending..." << endl;
		isjammed = -1;
	}
}

void rigidbody::rb_easy(double& phiH, double& phiL, int& check_rattlers, int &epconst, int nr, double dphi0, double Ktol, double &Utol, int t, bool &min) {
	// local variables
	int nbb, niso, pcsum, acsum;
	bool gr, oc, uc, ismin;
	double dphi = dphi0;

	// initialize bools
	gr = 0;
	oc = 0;
	uc = 0;

	// get number of contacts
	nbb = N - nr;
	niso = DOF * nbb - NDIM + 1;
	acsum = 0.5*get_ac_sum();
	pcsum = get_c_sum();

	// check bool conditions
	gr = (U < Utol);
	oc = (U > 2 * Utol && K < Ktol && epconst == 1 && acsum > 0 && !min);
	uc = (U < Utol && epconst == 1 && !min);
	ismin = (epconst == 1 && min);

	// increase/decrease phi based on U and K
	if (phiH < 0) {
		if (gr) {
			check_rattlers = 0;
			rb_scale(phi + dphi);
		}
		else if (oc && epconst == 1) {
			// phiH = phi;
			// dphi = -0.5*dphi0;
			// check_rattlers = 1;

			// cout << endl;
			// cout << "phiH 1st set at nt = " << t << endl;
			// monitor_scale(phi + dphi, phiL, phiH);

			// if NLCL, change update check (particles don't move, don't need to check as often)
			// if (NCL > -1)
			// 	nnupdate *= 50;

			// set min = 1, end
			min = true;
		}

		// Minimized, probably not as close as one would want but hey! points for trying!
		if (ismin){
			cout << endl;
			cout << endl;
			cout << "Found minimized state that's close enough!" << endl;
			cout << "Writing config at t = " << t*dt << endl;
			cout << "Writing config at nt = " << t << endl;
			cout << "Final U = " << U << endl;
			cout << "Final K = " << K << endl;
			cout << "Final phi = " << setprecision(6) << phi << endl;
			cout << "Final pcsum = " << pcsum << endl;
			cout << "Final acsum = " << acsum << endl;
			cout << "Final niso = " << niso << endl;
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

			cout << endl;
			cout << endl;
			isjammed = 1;
		}
	}
	else {
		if (phiL < 0) {

			// if still overcompressed, decrease again
			if (oc && epconst == 1) {
				phiH = phi;
				dphi = -0.25*dphi0;				
				cout << endl;
				cout << "still overcompressed..." << endl;
				cout << "phiH set at nt = " << t << endl;
				monitor_scale(phi + dphi, phiL, phiH);
			}

			// if undercompressed, increase packing fraction, done!
			if (uc) {
				phiL = phi;
				dphi = 0.5 * (phiH + phiL) - phi;

				cout << endl;
				cout << "relaxation found!" << endl;
				cout << "phiL set at nt = " << t << endl;
				monitor_scale(phi + dphi, phiL, phiH);
			}

		}
		else {

			// if undercompressed, root search up
			if (uc) {
				phiL = phi;
				dphi = 0.5 * (phiH + phiL) - phi;

				cout << endl;
				cout << "relaxed state found!" << endl;
				cout << "phiL set at nt = " << t << endl;
				monitor_scale(phi + dphi, phiL, phiH);
			}

			// if overcompressed, minimize and get out of there
			if (oc) {
				cout << endl;
				cout << "***************************" << endl;
				cout << "***************************" << endl;
				cout << "overcompressed state found!" << endl;
				cout << "min set to true at nt = " << t << endl;
				cout << "***************************" << endl;
				cout << "***************************" << endl;
				min = true;
				epconst = 0;
			}

			// STRICTER CONDITION
			// // Minimized, probably not as close as one would want but hey! points for trying!
			// if (ismin){
			// 	cout << endl;
			// 	cout << endl;
			// 	cout << "Found minimized state that's close enough!" << endl;
			// 	cout << "Writing config at t = " << t*dt << endl;
			// 	cout << "Writing config at nt = " << t << endl;
			// 	cout << "Final U = " << U << endl;
			// 	cout << "Final K = " << K << endl;
			// 	cout << "Final phi = " << setprecision(6) << phi << endl;
			// 	cout << "Final pcsum = " << pcsum << endl;
			// 	cout << "Final acsum = " << acsum << endl;
			// 	cout << "Final niso = " << niso << endl;
			// 	cout << "Final niso max = " << DOF*N - NDIM + 1 << endl;
			// 	cout << "Final rattler # = " << nr << endl;
			// 	cout << "Final contacts:" << endl;
			// 	cout << "pc:" << endl;
			// 	for (int i = 1; i < N + 1; i++) {
			// 		cout << setw(6) << pc[i - 1];
			// 		if (i % 10 == 0)
			// 			cout << endl;
			// 	}
			// 	cout << endl;
			// 	cout << "ac:" << endl;
			// 	for (int i = 1; i < N + 1; i++) {
			// 		cout << setw(6) << ac[i - 1];
			// 		if (i % 10 == 0)
			// 			cout << endl;
			// 	}
			// 	cout << endl;
			// 	if (N < 40) {
			// 		cout << "Contact Matrix:" << endl;
			// 		print_c_mat();
			// 	}
			// 	if (xyzobj.is_open())
			// 		rigidbody_xyz();

			// 	cout << endl;
			// 	cout << endl;
			// 	isjammed = 1;
			// }
			
		}
	}
}
















/*

==================================

			PRINTING

==================================

*/

// PRINT methods
void rigidbody::print_data() {
	// throw error if file not opened
	if (!configobj.is_open()) {
		cout << "ERROR: config file not opened!" << endl;
		throw "ERROR: config file not opened!";
	}

	// local variables
	int w, p, i, j, d;
	w = 14;
	p = 6;

	// print stat info
	configobj << setw(w) << "isjammed: " << isjammed << endl;
	configobj << setw(w) << "N: " << N << endl;
	configobj << setw(w) << "L: ";
	for (d = 0; d < NDIM; d++)
		configobj << setw(w) << L[d];
	configobj << endl;
	configobj << setw(w) << setprecision(16) << "phi: " << phi << endl;
	configobj << setw(w) << "U: " << U << endl;
	configobj << setw(w) << "K: " << K << endl;
	configobj << setw(w) << "pc sum: " << get_c_sum() << endl;
	configobj << setw(w) << "ac sum: " << 0.5 * (get_ac_sum()) << endl;
	configobj << setw(w) << "nr: " << nr << endl;
	configobj << setw(w) << "niso: " << DOF*(N - nr) - NDIM + 1 << endl;
	configobj << setw(w) << "niso max: " << DOF*N - NDIM + 1 << endl;
	configobj << setw(w) << "pc: ";
	print_pc(configobj, round(w / 2));
	configobj << setw(w) << "ac: ";
	print_ac(configobj, round(w / 2));
	configobj << "contact matrix: ";
	print_c_data(configobj);
	configobj << setw(w) << "configuration: " << endl;

	// update euler angles, given quaternions
	update_euler();

	// update print width
	w = 30;
	p = 16;

	// print header
	configobj << setw(w) << "id";
	configobj << setw(w) << "arad";
	for (d = 0; d < NDIM; d++)
		configobj << setw(w) << "ax" << d << "";
	configobj << setw(w) << "prad";
	for (d = 0; d < NDIM; d++)
		configobj << setw(w) << "px" << d << "";
	for (d = 0; d < NDIM; d++)
		configobj << setw(w) << "Inn" << d << "";
	configobj << setw(w) << "phi";
	configobj << setw(w) << "theta";
	configobj << setw(w) << "psi";
	configobj << setw(w) << "m";
	configobj << setw(w) << "Na";
	configobj << endl;
	for (i = 0; i < w * (3 * NDIM + 8); i++)
		configobj << "=";
	configobj << endl;

	// print data
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			configobj << setw(w) << i;
			configobj << setw(w) << setprecision(p) << ar[i][j];
			for (d = 0; d < NDIM; d++)
				configobj << setw(w) << setprecision(p) << xW[i][j][d]+x[i][d];
			configobj << setw(w) << setprecision(p) << r[i];
			for (d = 0; d < NDIM; d++)
				configobj << setw(w) << setprecision(p) << x[i][d];
			for (d = 0; d < NDIM; d++)
				configobj << setw(w) << setprecision(p) << Inn[i][d];
			configobj << setw(w) << setprecision(p) << eulang1[i];
			configobj << setw(w) << setprecision(p)	<< eulang2[i];
			configobj << setw(w) << setprecision(p) << eulang3[i];
			configobj << setw(w) << setprecision(p) << m[i];
			configobj << setw(w) << setprecision(p) << Na[i];
			configobj << endl;
		}
	}
}

void rigidbody::print_stat() {
	// throw error if file not opened
	if (!statobj.is_open()) {
		cout << "ERROR: config file not opened!" << endl;
		throw "ERROR: config file not opened!";
	}

	// local variables
	int w, p, i, d;
	w = 14;
	p = 6;

	// print stat info
	statobj << setw(w) << "isjammed: " << isjammed << endl;
	statobj << setw(w) << "N: " << N << endl;
	statobj << setw(w) << "L: ";
	for (d = 0; d < NDIM; d++)
		statobj << setw(w) << L[d];
	statobj << endl;
	statobj << setw(w) << "phi: " << phi << endl;
	statobj << setw(w) << "U: " << U << endl;
	statobj << setw(w) << "K: " << K << endl;
	statobj << setw(w) << "pc sum: " << get_c_sum() << endl;
	statobj << setw(w) << "ac sum: " << 0.5 * (get_ac_sum()) << endl;
	statobj << setw(w) << "nr: " << nr << endl;
	statobj << setw(w) << "niso: " << DOF*(N - nr) - NDIM + 1 << endl;
	statobj << setw(w) << "niso max: " << DOF*N - NDIM + 1 << endl;
	statobj << setw(w) << "pc: ";
	print_pc(statobj, round(w / 2));
	statobj << setw(w) << "ac: ";
	print_ac(statobj, round(w / 2));
	statobj << "contact matrix: ";
	print_c_data(statobj);
	statobj << endl;
}

void rigidbody::print_ac(ofstream& obj, int w) {
	int i;
	for (i = 0; i < N; i++)
		obj << setw(w) << ac[i];
	obj << endl;
}

void rigidbody::print_config() {
	// throw error if file not opened
	if (!configobj.is_open()) {
		cout << "ERROR: config file not opened!" << endl;
		throw "ERROR: config file not opened!";
	}

	// local variables
	int w, p, i, j, d;
	w = 30;
	p = 16;

	// update euler angles, given quaternions
	update_euler();

	// print basic info
	configobj << setw(w) << N << endl;
	for (d = 0; d < NDIM; d++)
		configobj << setw(w) << L[d];
	configobj << endl;
	configobj << setw(w) << setprecision(p) << phi << endl;

	// print header
	configobj << setw(w) << "id";
	configobj << setw(w) << "arad";
	for (d = 0; d < NDIM; d++)
		configobj << setw(w) << "ax" << d << "";
	configobj << setw(w) << "prad";
	for (d = 0; d < NDIM; d++)
		configobj << setw(w) << "px" << d << "";
	for (d = 0; d < NDIM; d++)
		configobj << setw(w) << "Inn" << d << "";
	configobj << setw(w) << "phi";
	configobj << setw(w) << "theta";
	configobj << setw(w) << "psi";
	configobj << setw(w) << "m";
	configobj << setw(w) << "Na";
	configobj << endl;
	for (i = 0; i < w * (3 * NDIM + 8); i++)
		configobj << "=";
	configobj << endl;

	// print data
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			configobj << setw(w) << i;
			configobj << setw(w) << setprecision(p) << ar[i][j];
			for (d = 0; d < NDIM; d++)
				configobj << setw(w) << setprecision(p) << xW[i][j][d]+x[i][d];
			configobj << setw(w) << setprecision(p) << r[i];
			for (d = 0; d < NDIM; d++)
				configobj << setw(w) << setprecision(p) << x[i][d];
			for (d = 0; d < NDIM; d++)
				configobj << setw(w) << setprecision(p) << Inn[i][d];
			configobj << setw(w) << setprecision(p) << eulang1[i];
			configobj << setw(w) << setprecision(p)	<< eulang2[i];
			configobj << setw(w) << setprecision(p) << eulang3[i];
			configobj << setw(w) << setprecision(p) << m[i];
			configobj << setw(w) << setprecision(p) << Na[i];
			configobj << endl;
		}
	}
}

void rigidbody::rigidbody_md_monitor() {
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
	cout << "sum c = " << get_c_sum() << endl;
	cout << "sum ac = " << 0.5 * (get_ac_sum()) << endl;
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
		enobj << setw(12) << K/N;
		enobj << setw(12) << Krot/N;
		enobj << setw(12) << U/N + K/N;
		enobj << setw(12) << alpha;
		enobj << endl;
	}	

	// output to xyz file if open
	if (xyzobj.is_open()) {
		cout << "Printing XYZ" << endl;
		rigidbody_xyz();
		cout << endl;
	}
}

void rigidbody::monitor_scale(double phinew, double phiL, double phiH) {
	cout << "phiH = " << phiH << endl;
	cout << "phiL = " << phiL << endl;
	cout << "dphi = " << phinew - phi << endl;
	cout << "old phi = " << phi << endl;
	rb_scale(phinew);
	cout << "new phi = " << phi << endl;
	cout << "U = " << U/N << endl;
	cout << "K = " << K/N << endl;
	cout << "c = " << 0.5*get_ac_sum() << endl;
	cout << "conitinuing root search..." << endl;
	cout << endl;
}

void rigidbody::rigidbody_xyz() {
	int i, j, k, nf, d, dd, w, n_found;
	double rmin, rtmp;
	w = 20;

	// get min radius (for scaling)
	rmin = 1e20;
	for (j = 0; j < Na[0]; j++) {
		if (ar[0][j] < rmin)
			rmin = ar[0][j];
	}

	// print xyz header
	xyzobj << get_Natot() << endl;
	xyzobj << "Lattice=\"";
	for (d = 0; d < NDIM; d++) {
		for (dd = 0; dd < NDIM; dd++) {
			if (dd == d)
				xyzobj << L[d];
			else
				xyzobj << " 0.0 ";
		}
	}
	xyzobj << "\" ";
	xyzobj << '\t';
	xyzobj << "Properties=particleID:R:1:MoleculeType:S:1:species:S:1:pos:R:3:radius:R:1:Velocity:R:3:AngularVelocity:R:3" << endl;

	// rotate rotational velocity to world frame
	Quaternion q1;
	Quaternion q2;
	Quaternion wtmp;
	for (i=0; i<N; i++){
		// rotate rel M frame to rel W frame
		wtmp.set_s(0.0);
		wtmp.set_x(wM[i][0]);
		wtmp.set_y(wM[i][1]);
		wtmp.set_z(wM[i][2]);

		q1 *= q[i];
		q2 = wtmp % q1;
		wtmp = q[i] % q2;

		wW[i][0] = wtmp.get_x();
		wW[i][1] = wtmp.get_y();
		wW[i][2] = wtmp.get_z();
	}

	// print body
	n_found = 0;
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			// print particle (i.e. residue) id
			xyzobj << setw(w) << i;

			// print res name
			switch (Na[i]){
				case 10: xyzobj << setw(w) << "A"; break;
				case 16: xyzobj << setw(w) << "V"; break;
				case 17: xyzobj << setw(w) << "M"; break;
				case 19: xyzobj << setw(w) << "L"; break;
				case 20: xyzobj << setw(w) << "P"; break;
				default: xyzobj << setw(w) << "X"; break;
			}

			// determine atomic type based on radius
			rtmp = ar[i][j] / rmin;
			if (rtmp < 1.12) {
				xyzobj << setw(w) << 'H';
			}			
			else if (rtmp < 1.36) {
				if (n_found == 0) {
					xyzobj << setw(w) << 'N' << i;
					n_found = 1;
				}
				else {
					xyzobj << setw(w) << 'C' << i;
					n_found = 0;
				}
			}
			else if (rtmp < 1.44)
				xyzobj << setw(w) << 'O';
			else if (rtmp < 1.56)
				xyzobj << setw(w) << 'C';
			else
				xyzobj << setw(w) << 'S';				

			// print com positions
			xyzobj << setw(w) << xW[i][j][0] + x[i][0];
			xyzobj << setw(w) << xW[i][j][1] + x[i][1];
			xyzobj << setw(w) << xW[i][j][2] + x[i][2];

			// print atomic radii
			xyzobj << setw(w) << ar[i][j];

			// print com translational velocities
			xyzobj << setw(w) << v[i][0];
			xyzobj << setw(w) << v[i][1];
			xyzobj << setw(w) << v[i][2];

			// print rigid body angular velocity in world frame
			xyzobj << setw(w) << wW[i][0];
			xyzobj << setw(w) << wW[i][1];
			xyzobj << setw(w) << wW[i][2];

			// print new line
			xyzobj << endl;
		}
	}
}

void rigidbody::rigidbody_xyz(int p1, int p2) {
	int i, j, k, nf, d, dd, w, n_found;
	double rmin, rtmp;
	w = 20;

	// get min radius (for scaling)
	rmin = 1e20;
	for (j = 0; j < Na[0]; j++) {
		if (ar[0][j] < rmin)
			rmin = ar[0][j];
	}

	// print xyz header
	xyzobj << get_Natot() << endl;
	xyzobj << "Lattice=\"";
	for (d = 0; d < NDIM; d++) {
		for (dd = 0; dd < NDIM; dd++) {
			if (dd == d)
				xyzobj << L[d];
			else
				xyzobj << " 0.0 ";
		}
	}
	xyzobj << "\" ";
	xyzobj << '\t';
	xyzobj << "Properties=species:S:1:pos:R:" <<  NDIM << ":radius:R:1" << endl;

	// print body
	n_found = 0;
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			// determine atomic type
			rtmp = ar[i][j] / rmin;
			if (i==p1)
				xyzobj << setw(w) << 'X';
			else if (i==p2)
				xyzobj << setw(w) << 'Y';
			else{
				if (rtmp < 1.12)					
					xyzobj << setw(w) << 'H';
				else if (rtmp < 1.36) {
					if (n_found == 0) {
						xyzobj << setw(w) << 'N';
						n_found = 1;
					}
					else {
						xyzobj << setw(w) << 'C';
						n_found = 0;
					}
				}
				else if (rtmp < 1.44)
					xyzobj << setw(w) << 'O';
				else if (rtmp < 1.56)
					xyzobj << setw(w) << 'C';
				else
					xyzobj << setw(w) << 'S';				
			}
			

			xyzobj << setw(w) << xW[i][j][0] + x[i][0];
			xyzobj << setw(w) << xW[i][j][1] + x[i][1];
			xyzobj << setw(w) << xW[i][j][2] + x[i][2];
			xyzobj << setw(w) << ar[i][j];
			xyzobj << endl;
		}
	}
}

void rigidbody::rigidbody_print_vars() {
	int i, j, d;

	cout << "** Printing packing vars:" << endl;
	print_vars();

	cout << endl << "** rigid body arrays:" << endl << endl;
	cout << "xW: " << endl;
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			for (d = 0; d < NDIM; d++)
				cout << setw(20) << setprecision(6) << "xW[" << i << "][" << j << "][" << d << "] = " << xW[i][j][d];
			cout << endl;
		}
	}
	cout << endl;
	cout << "xM: " << endl;
	for (i = 0; i < N; i++) {
		for (j = 0; j < Na[i]; j++) {
			for (d = 0; d < NDIM; d++)
				cout << setw(20) << setprecision(6) << "xM[" << i << "][" << j << "][" << d << "] = " << xM[i][j][d];
			cout << endl;
		}
	}
	cout << endl;
	cout << "q: " << endl;
	for (i = 0; i < N; i++) {
		cout << "q[" << i << "]" << endl;
		q[i].print_vals();
		cout << endl;
	}
	cout << endl;
	cout << "euler angles: " << endl;
	for (i = 0; i < N; i++) {
		cout << setw(20) << setprecision(6) << "alpha[" << i << "] = " << eulang1[i];
		cout << setw(20) << setprecision(6) << "beta[" << i << "] = " << eulang2[i];
		cout << setw(20) << setprecision(6) << "gamma[" << i << "] = " << eulang3[i];
		cout << endl;
	}
	cout << endl;
}














