// quaternion test

#include "Quaternion.h"
#include <iostream>

int main(){
	double q1s,q1x,q1y,q1z;
	double q2s,q2x,q2y,q2z;

	// initialize q1
	q1s = -1.0;
	q1x = 1.0;
	q1y = 2.0;
	q1z = 3.0;
	Quaternion q1(q1s,q1x,q1y,q1z);

	// initialize q2
	q2s = 3.0;
	q2x = 1.0;
	q2y = 2.0;
	q2z = -1.0;
	Quaternion q2(q2s,q2x,q2y,q2z);

	// do some operations, print vals
	cout << "q1: " << endl;
	q1.print_vals();
	cout << endl;

	cout << "q2: " << endl;
	q2.print_vals();
	cout << endl;

	// make q3 from addition of q1 and q2
	Quaternion q3;
	q3 = q1+q2;
	
	cout << "q3 = q1+q2: " << endl;
	q3.print_vals();
	cout << endl;

	// dot product	
	double val;
	val = q1 * q2;
	cout << "dot product = " << val << endl << endl;

	// make q3 from quaternion multiplication of q1 and q2
	q3 = q1 % q2;

	cout << "q3 = q1 % q2: " << endl;
	q3.print_vals();
	cout << endl;

	// q4: vector in quaternion form
	Quaternion q4(0,1,3,2);
	q1.normalize();
	q2 *= q1;
	q3 = q2 % (q4 % q1);
	cout << "q3 = q1* q4 q1:" << endl;
	q3.print_vals();
	cout << endl;



	
	return 0;
}