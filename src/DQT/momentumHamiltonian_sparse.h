#ifndef SPINCHAINED_MOMENTUMHAMILTONIAN_SPARSE_H
#define SPINCHAINED_MOMENTUMHAMILTONIAN_SPARSE_H

#include <vector>
#include <list>
#include <complex>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "../ED/helperFunctions.h"
#include "../ED/momentumHamiltonian.h"

Eigen::SparseMatrix<std::complex<double>> momentumHamiltonian_sparse(double J_ratio, int N);

Eigen::SparseMatrix<std::complex<double>> spinOp2_momentum_sparse(int N);

Eigen::SparseMatrix<std::complex<double>> spinOp2_momentum_sparse_large_N(int N);

Eigen::SparseMatrix<std::complex<double>> spinOpS2_momentum_sparse_m(int N, int n_up);

std::vector<Eigen::SparseMatrix<std::complex<double>>> momentumHamiltonian_sparse_blocks(double J_ratio, int N, int m_lower, int m_upper);


#endif //SPINCHAINED_MOMENTUMHAMILTONIAN_SPARSE_H
