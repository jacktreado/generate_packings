/*

	Methods implementation
	for trajectorty analysis of 
	rigidbody class

	BY Jack Treado

*/

#include "Quaternion.h"
#include "rigidbody.h"
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

const double PI = atan(1)*4;


void rigidbody::free_md(double tmp0, int NT, int nnu) {
	int t;

	// initialize velocities
	this->rand_vel_init(tmp0);

	// if energy output open, output dt
	enobj << dt << endl;

	// update nearest neighbor update
	nnupdate = nnu;

	cout << "===== BEGINNING TIME LOOP ======" << endl;
	cout << "NT = " << NT << endl;
	cout << "dt = " << dt << endl;

	for (t = 0; t < NT; t++) {
		// if NLCL, update
		if (NCL > -1 && t % nnupdate == 0) {
			this->update_cell();
			this->update_neighborlist();
		}

		// advance quaternions, positions
		this->verlet_first();

		// update forces
		this->force_update();

		// advance angular momentum
		this->verlet_second();

		// output information
		if (t % plotskip == 0) {
			this->monitor_header(t);
			this->rigidbody_md_monitor();
		}
	}
}