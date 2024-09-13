#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main(){
    const int nbos = 2;
    const int nsites = 2;
    const double t = -1.0;
    const double omega = 1.0;
    const double g = 0.5;
    
    MatrixXd b = MatrixXd::Zero(nbos, nbos);

    for (int j = 1; j<nbos; ++j){
        b(j-1,j) = sqrt(j);
    }

    //cout<<"Matrix b: \n" <<b <<endl;
    
    MatrixXd b_dagger = b.transpose();
    MatrixXd sum_b_b_dagger = b_dagger + b;

    MatrixXd hop = MatrixXd::Zero(nsites,nsites);
    for (int j = 1; j < nsites; ++j) {
	    hop(j-1,j) = t;
	    hop(j,j-1) = t;
    }

    cout << "HOPPING\n" << hop <<endl; 
    return 0;
}
