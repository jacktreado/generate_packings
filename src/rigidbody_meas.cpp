/*

	Measure methods implementation
	for rigidbody class

	BY Jack Treado
	10/01/2018

*/

#include "rigidbody.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <vector>




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