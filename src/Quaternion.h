/* * * * * * * * * * * * * * * * *

	Quaternion Class
	by Jack Treado, 08/24/17

	Store scalar and vector quantities in Quaternion objects

 * * * * * * * * * * * * * * * * */

#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>
using namespace std;

class Quaternion
{
	private:	
		double s; 		// scalar value 
		double x;
		double y;
		double z;
	public:

		// constructors
		Quaternion();
		Quaternion(double s, double x, double y, double z);
		Quaternion(const Quaternion& qinit);

		// destructor
		~Quaternion();

		// setters
		void set_s(double s2){ s = s2; };
		void set_x(double x2){ x = x2; };
		void set_y(double y2){ y = y2; };
		void set_z(double z2){ z = z2; };
		void set_rand(int seed);

		// getters	
		double get_norm() const;
		double get_s() const { return s; };
		double get_x() const { return x; };
		double get_y() const { return y; };
		double get_z() const { return z; };

		// overloaded operators
		void operator= (const Quaternion &qright);							// assigment operator
		void operator*= (const Quaternion &qright);							// assign conjugate of rhs
		void operator-= (const Quaternion &qright);							// assign inverse of rhs
		Quaternion operator+ (const Quaternion &qright);					// addition operator
		Quaternion operator- (const Quaternion &qright);					// subtraction operator
		Quaternion operator* (const double dright);							// product with scalar (double)
		Quaternion operator* (const int iright);							// product with scalar (integer)
		double operator* (const Quaternion &qright);						// dot product operator
		Quaternion operator% (const Quaternion &qright);					// vector product operator

		// other operations		
		void conjugate();
		void invert();
		void normalize();

		// Printers & test
		void print_vals();
};


#endif