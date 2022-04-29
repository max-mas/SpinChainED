//
// Created by mmaschke on 4/25/22.
//

#ifndef HEISENBERG_CHAIN_1D_C_MAIN_H
#define HEISENBERG_CHAIN_1D_C_MAIN_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <list>

#include <cmath>
#include <complex>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "helperFunctions.h"

std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXcd> & H_list, int n);

std::vector<double> getEnergiesFromBlocks(const std::list<std::list<Eigen::MatrixXcd>> & H_list, int N);

Eigen::MatrixXcd blkdiag(const std::vector<Eigen::MatrixXcd>& matrix_list, int totalSize);

Eigen::MatrixXd blkdiag(const std::vector<Eigen::MatrixXd> & matrix_list, int totalSize);



#endif //HEISENBERG_CHAIN_1D_C_MAIN_H
