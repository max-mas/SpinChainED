
#include "matrixExpIteration.h"

using Eigen::SparseMatrix;
using Eigen::VectorXcd;
using Eigen::VectorXd;

using std::complex;

void iterateState_beta(const SparseMatrix<complex<double>> & H, VectorXcd & psi, double dBeta) {
    VectorXcd k1 = -0.5 * H *  psi;
    VectorXcd k2 = -0.5 * H * (psi + k1 * 0.5 * dBeta);
    VectorXcd k3 = -0.5 * H * (psi + k2 * 0.5 * dBeta);
    VectorXcd k4 = -0.5 * H * (psi + k3 * dBeta);

    psi += dBeta/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
    //psi.normalize();
}


