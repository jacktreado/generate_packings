/*

	Test file to check use of
	Spectra and Eigen eigenvalue
	solvers

*/

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <iostream>

using namespace std;


int main(){

	// We are going to solve the generalized eigenvalue problem A * x = lambda * B * x
    const int n = 10;


    // Define the A matrix
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd A = M + M.transpose();

    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    Eigen::SparseMatrix<double> B(n, n);
    B.reserve(Eigen::VectorXi::Constant(n, 3));
    for(int i = 0; i < n; i++)
    {
        B.insert(i, i) = 2.0;
        if(i > 0)
            B.insert(i - 1, i) = 1.0;
        if(i < n - 1)
            B.insert(i + 1, i) = 1.0;
    }    

    // Verify results using the generalized eigen solver in Eigen
    Eigen::MatrixXd Bdense = B;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Bdense);
    cout << "Generalized eigenvalues (verify):\n" << es.eigenvalues() << endl;
    cout << "Generalized eigenvectors (verify):\n" << es.eigenvectors().rightCols(n).topRows(n) << endl;


	return 0;
}