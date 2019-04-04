// CLASS AND METHODS DEFINTION FOR vec3 CLASS

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#ifndef VEC3_H
#define VEC3_H

class vec3{
private:

	// vector components
	double X;
	double Y;
	double Z;

public:

	// default constructor
	vec3(){ X=0.0; Y=0.0; Z=0.0; }

	// overloaded constructors
	vec3(double s) { X=s; Y=s; Z=s; };
	vec3(double x, double y, double z) { X=x; Y=y; Z=z; };
	vec3(std::vector<double>& v){
		// check input size
		if (v.size() != 3){
			std::cout << "ERROR: input vector wrong size" << std::endl;
			throw;
		}

		// if vector is ok size, assign values to private variables
		X = v.at(0);
		Y = v.at(1);
		Z = v.at(2);
	};
	vec3(double arr[]){
		X = arr[0];
		Y = arr[1];
		Z = arr[2];
	}

	// setters
	void set_x (double x) { X=x; };
	void set_y (double y) { Y=y; };
	void set_z (double z) { Z=z; };
	void reset (double arr[]) {
		X = arr[0];
		Y = arr[1];
		Z = arr[2];
	};

	// getters
	double at(int i) const;
	double get_x() const { return X; };
	double get_y() const { return Y; };
	double get_z() const { return Z; };
	double get_norm() { 
		double s = 0.0; 
		s = sqrt(X*X + Y*Y + Z*Z); 
		return s; 
	};

	// self-operations
	void normalize(){
		double s = 0.0;
		s = this->get_norm();

		X /= s;
		Y /= s;
		Z /= s;
	}

	// overloaded operators
	void operator= (const vec3 &vright);			// assigment operator
	vec3 operator+ (const vec3 &vright);			// addition operator
	vec3 operator- (const vec3 &vright);			// subtraction operator
	vec3 operator* (const double dright);			// product with scalar (double)
	vec3 operator* (const int iright);				// product with scalar (integer)
	double operator* (const vec3 &vright);			// dot product operator
	vec3 operator% (const vec3 &vright);			// vector product operator

	// print methods
	void print(){
		std::cout << std::setw(10) << X;
		std::cout << std::setw(10) << Y;
		std::cout << std::setw(10) << Z;
		std::cout << std::endl;
	};

};

#endif










