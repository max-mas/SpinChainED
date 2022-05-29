
#ifndef HEISENBERG_CHAIN_1D_C__NAIVEHAMILTONIAN_SPARSE_H
#define HEISENBERG_CHAIN_1D_C__NAIVEHAMILTONIAN_SPARSE_H

#include <list>
#include <vector>

#include <Eigen/Sparse>

#include "../ED/helperFunctions.h"

Eigen::SparseMatrix<double> naiveHamiltonian_sparse(double J_ratio, int N);


#endif //HEISENBERG_CHAIN_1D_C__NAIVEHAMILTONIAN_SPARSE_H
