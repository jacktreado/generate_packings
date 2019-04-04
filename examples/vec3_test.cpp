// test main cpp file for vec3 class

#include <iostream>
#include <iomanip>
#include "vec3.h"

#define L 3

using namespace std;

int main(){
	// local variables
	vec3 v1(1,0,0);
	vec3 v2(0,1,0);

	// check methods
	cout << endl;
	cout << "vector 1: ";
	v1.print();
	cout << endl;

	cout << endl;
	cout << "vector 2: ";
	v2.print();
	cout << endl;


	// do vector dot and cross product;
	cout << endl;
	cout << "v1 * v2 = ";
	vec3 v3;
	v3 = v1 * v2;
	v3.print();
	cout << endl;

	cout << endl;
	cout << "v1 x v2 = ";
	v3 = v1 % v2;
	v3.print();
	cout << endl;

	// change v3
	v3.set_y(3.0);

	// output
	cout << endl;
	cout << "vector 3 now = ";
	v3.print();
	cout << endl;

	// normalize
	cout << endl;
	v3.normalize();
	cout << "after normalizing, vector 3 = ";
	v3.print();

	// to check, print norm
	cout << "new norm of v3 = " << v3.get_norm();
	cout << endl;



	// test use of array input
	double testArr[L][L];
	double testColArr[L];

	cout << endl << endl << "testing array passing:" << endl;
	for (int i=0; i<L; i++){
		for (int j=0; j<L; j++){
			testArr[i][j] = i-2*j;
			cout << setw(6) << testArr[i][j];
		}
		cout << endl;

		// fill up column array
		testColArr[i] = 100*i - 99;
	}

	cout << "passing array as two vectors: " << endl;
	vec3 varr1(testArr[0]);
	vec3 varr2(testArr[1]);
	vec3 vColArr(testColArr);

	cout << "printing..." << endl;
	varr1.print();
	varr2.print();
	vColArr.print();


	// end program
	return 0;
}