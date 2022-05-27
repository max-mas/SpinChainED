
#include "matrixExpIteration.h"

using Eigen::SparseMatrix;
using Eigen::VectorXcd;

using std::complex;

void iterateState_beta(const SparseMatrix<complex<double>> & H, VectorXcd & psi, double dBeta) {
    VectorXcd k1 = -complex<double>(0, 0.5) * H * psi;
    VectorXcd k2 = -complex<double>(0, 0.5) * H * (psi + 0.5 * dBeta * k1);
    VectorXcd k3 = -complex<double>(0, 0.5) * H * (psi + 0.5 * dBeta * k2);
    VectorXcd k4 = -complex<double>(0, 0.5) * H * (psi + dBeta * k3);

    psi += dBeta/6.0 * ( k1 + 2.0*k2 + 2.0*k3 + k4 );
    //psi.normalize();
}


