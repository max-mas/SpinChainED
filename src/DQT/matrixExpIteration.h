
#ifndef SPINCHAINED_MATRIXEXPITERATION_H
#define SPINCHAINED_MATRIXEXPITERATION_H

#include <complex>

#include <Eigen/Sparse>
#include <Eigen/Dense>

void iterateState_beta(const Eigen::SparseMatrix<std::complex<double>> & H, Eigen::VectorXcd & psi, double dBeta);


#endif //SPINCHAINED_MATRIXEXPITERATION_H
