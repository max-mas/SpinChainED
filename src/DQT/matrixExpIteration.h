/**
 * This file contains the beta iteration of the imaginary time Schrödinger equation using a 4th order Runge-Kutta
 * scheme.
 */

#ifndef SPINCHAINED_MATRIXEXPITERATION_H
#define SPINCHAINED_MATRIXEXPITERATION_H

#include <complex>

#include <Eigen/Sparse>
#include <Eigen/Dense>

// Implements RK4 for the imaginary time SE.
void iterateState_beta(const Eigen::SparseMatrix<std::complex<double>> & H, Eigen::VectorXcd & psi, double dBeta);

void iterateState_beta(const Eigen::SparseMatrix<double> & H, Eigen::VectorXcd & psi, double dBeta);

#endif //SPINCHAINED_MATRIXEXPITERATION_H
