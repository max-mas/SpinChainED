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

Eigen::MatrixXd naiveHamiltonian(double J_ratio, int N);

std::list<Eigen::MatrixXd> magnetizationHamiltonian(double J_ratio, int N);

std::list<std::list<Eigen::MatrixXcd>> momentumHamiltonian(double J_ratio, int N);

std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXcd> & H_list, int N);

std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXd> & H_list, int N);

std::vector<double> getEnergiesFromBlocks(const std::list<std::list<Eigen::MatrixXcd>> & H_list, int N);

Eigen::MatrixXcd blkdiag(const std::list<Eigen::MatrixXcd>& matrix_list, int totalSize);

std::list<Eigen::MatrixXcd> blkdiag(const std::list<std::list<Eigen::MatrixXcd>> & matrix_doubleList);

#endif //HEISENBERG_CHAIN_1D_C_MAIN_H
