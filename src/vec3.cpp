// implementation for vec3 class

#include "vec3.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// retrival
double vec3::at(int i) const{
	if (i > 2){
		std::cout << "trying to access past NDIM = 3, throwing error" << std::endl;
		std::cout << std::endl;
		throw;
	}
	else{
		if (i==0)
			return X;
		else if (i==1)
			return Y;
		else if (i==2)
			return Z;
		else
			return -1;
	}

}


// METHODS DEFINITION

// assignment
void vec3::operator= (const vec3 &vright){
	X = vright.get_x();
	Y = vright.get_y();
	Z = vright.get_z();
}

// addition
vec3 vec3::operator+ (const vec3 &vright){
	// locally initialize output object
	vec3 vecOut;

	// set each element in output object
	vecOut.set_x(X + vright.get_x());
	vecOut.set_y(Y + vright.get_y());
	vecOut.set_z(Z + vright.get_z());

	// return object
	return vecOut;
}

// subtraction
vec3 vec3::operator- (const vec3 &vright){
	// locally initialize output object
	vec3 vecOut;

	// set each element in output object
	vecOut.set_x(X - vright.get_x());
	vecOut.set_y(Y - vright.get_y());
	vecOut.set_z(Z - vright.get_z());

	// return object
	return vecOut;
}

// scalar (double) multiplication
vec3 vec3::operator* (const double dright){
	// locally initialize output object
	vec3 vecOut;

	// set each element in output object
	vecOut.set_x(X*dright);
	vecOut.set_y(Y*dright);
	vecOut.set_z(Z*dright);

	// return object
	return vecOut;
}

// scalar (integer) multiplication
vec3 vec3::operator* (const int iright){
	// locally initialize output object
	vec3 vecOut;

	// set each element in output object
	vecOut.set_x(X*iright);
	vecOut.set_y(Y*iright);
	vecOut.set_z(Z*iright);

	// return object
	return vecOut;
}


// scalar dot product
double vec3::operator* (const vec3 &vright){
	// local variable to be returned
	int val = 0.0;

	// set each element in output object
	val = X*vright.get_x();
	val += Y*vright.get_y();
	val += Z*vright.get_z();

	// return object
	return val;
}


// vector cross product
vec3 vec3::operator% (const vec3 &vright){
	// locally initialize output object
	vec3 vecOut;
	double valx, valy, valz;
	double vrx, vry, vrz;

	// locally store values in input object
	vrx = vright.get_x();
	vry = vright.get_y();
	vrz = vright.get_z();

	// get values of each element
	valx = Y*vrz - Z*vry;
	valy = Z*vrx - X*vrz;
	valz = X*vry - Y*vrx;

	// set values to each element in vector
	vecOut.set_x(valx);
	vecOut.set_y(valy);
	vecOut.set_z(valz);

	// return object
	return vecOut;
}

