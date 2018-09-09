/* * * * * * * * * * * * * * * * *

Quaternion Class Methods Implementation
by Jack Treado, 08/24/17

Store scalar and vector quantities in Quaternion objects


 * * * * * * * * * * * * * * * * */


#include "Quaternion.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <cmath>
using namespace std;


// CONSTRUCTORS & DESTRUCTORS
Quaternion::Quaternion(){
	s = 0;
	x = 0;
	y = 0;
	z = 0;
}

Quaternion::Quaternion(double s2, double x2, double y2, double z2){
	s = s2;
	x = x2;
	y = y2;
	z = z2;
}

Quaternion::Quaternion(const Quaternion& qinit){
	s = qinit.get_s();
	x = qinit.get_x();
	y = qinit.get_y();
	z = qinit.get_z();
}

Quaternion::~Quaternion(){}




// SETTERS AND GETTERS

// Set to random numbers
void Quaternion::set_rand(int seed){
	srand48(seed);
	s = drand48();
	x = drand48();
	y = drand48();
	z = drand48();
}

double Quaternion::get_norm() const{
	double val = sqrt(s*s + x*x + y*y + z*z);
	return val;
}




// OVERLOADED OPERATORs

// assignment
void Quaternion::operator= (const Quaternion &qright){
	s = qright.get_s();
	x = qright.get_x();
	y = qright.get_y();
	z = qright.get_z();
}

// conjugate assigment
void Quaternion::operator*= (const Quaternion &qright){
	s = qright.get_s();
	x = -qright.get_x();
	y = -qright.get_y();
	z = -qright.get_z();
}

// inversion assigment
void Quaternion::operator-= (const Quaternion &qright){
	double val;

	// get denominator (complex product of orig * conj)
	val = pow(qright.get_s(),2) - pow(qright.get_x(),2) - pow(qright.get_y(),2) - pow(qright.get_z(),2);

	// assign to lhs
	s = qright.get_s()/val;
	x = qright.get_x()/val;
	y = qright.get_y()/val;
	z = qright.get_z()/val;	
}

// addition
Quaternion Quaternion::operator+ (const Quaternion &qright){
	// declare object to be copied and returned
	Quaternion qout;

	// add values to this
	qout.set_s(s+qright.get_s());
	qout.set_x(x+qright.get_x());
	qout.set_y(y+qright.get_y());
	qout.set_z(z+qright.get_z());
	return qout;
}

// subtraction
Quaternion Quaternion::operator- (const Quaternion &qright){
	// declare object to be copied and returned
	Quaternion qout;

	// add values to this
	qout.set_s(s-qright.get_s());
	qout.set_x(x-qright.get_x());
	qout.set_y(y-qright.get_y());
	qout.set_z(z-qright.get_z());
	return qout;
}

// overloaded multiplication/division with an int
Quaternion Quaternion::operator* (const double dright){
	// declare output object
	Quaternion qout;

	// add values to this
	qout.set_s(s*dright);
	qout.set_x(x*dright);
	qout.set_y(y*dright);
	qout.set_z(z*dright);
	return qout;
}

// overloaded multiplication/division with an int
Quaternion Quaternion::operator* (const int iright){
	// declare output object
	Quaternion qout;

	// add values to this
	qout.set_s(s*iright);
	qout.set_x(x*iright);
	qout.set_y(y*iright);
	qout.set_z(z*iright);
	return qout;
}

// dot product
double Quaternion::operator* (const Quaternion &qright){
	double dotp;

	dotp = s*qright.get_s() + x*qright.get_x() + y*qright.get_y() + z*qright.get_z();
	return dotp;
}

// quaternion-quaternion multiplication
Quaternion Quaternion::operator% (const Quaternion &qright){
	Quaternion qout;
	double vecps,vecpx,vecpy,vecpz;
	double crosspx,crosspy,crosspz;
	double q1s,q1x,q1y,q1z;
	double q2s,q2x,q2y,q2z;

	// get s,x,y,z for qleft
	q1s = s;
	q1x = x;
	q1y = y;
	q1z = z;	

	// get s,x,y,z for qright
	q2s = qright.get_s();
	q2x = qright.get_x();
	q2y = qright.get_y();
	q2z = qright.get_z();

	// get cross product components
	crosspx = q1y*q2z - q1z*q2y;
	crosspy = -(q1x*q2z - q1z*q2x);
	crosspz = q1x*q2y - q1y*q2x;

	vecps = q1s*q2s - (q1x*q2x + q1y*q2y + q1z*q2z);
	vecpx = q1s*q2x + q2s*q1x + crosspx;
	vecpy = q1s*q2y + q2s*q1y + crosspy;
	vecpz = q1s*q2z + q2s*q1z + crosspz;

	// set values of output quaternion
	qout.set_s(vecps);
	qout.set_x(vecpx);
	qout.set_y(vecpy);
	qout.set_z(vecpz);

	return qout;
}




// FUNCTIONS

void Quaternion::conjugate(){
	x *= -1;
	y *= -1;
	z *= -1;
}

void Quaternion::invert(){
	double val;
	val = s*s - x*x - y*y - z*z;

	// divide by conjugate product
	s /= val;
	x /= val;
	y /= val;
	z /= val;
}

void Quaternion::normalize(){
	double val = this->get_norm();

	// divide by norm
	s /= val;
	x /= val;
	y /= val;
	z /= val;
}


// Printer
void Quaternion::print_vals(){
	cout << "printing quaternions values: " << endl;
	cout << "s = " << s << endl;
	cout << "x = " << x << endl;
	cout << "y = " << y << endl;
	cout << "z = " << z << endl;
}








