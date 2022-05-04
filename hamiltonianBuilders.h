//
// Created by mmaschke on 03/05/22.
//

#ifndef SPINCHAINED_HAMILTONIANBUILDERS_H
#define SPINCHAINED_HAMILTONIANBUILDERS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <list>

#include <stdexcept>
#include <fstream>
#include <string>

#include <cmath>
#include <complex>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "helperFunctions.h"

Eigen::MatrixXd spinOperator_sq(int N);

Eigen::MatrixXd naiveHamiltonian(double J_ratio, int N);

void setHElement_naive(double J_ratio, int N, Eigen::MatrixXd &H, int a);

std::list<Eigen::MatrixXd> magnetizationHamiltonian(double J_ratio, int N);

std::vector<int> getStates_m(int N, int n_up);

void setHElement_magnetization(double J_ratio, int N, std::list<Eigen::MatrixXd> &H_subspace_list,
                               const std::vector<int> &s_vector_m, int k);

std::list<std::list<Eigen::MatrixXcd>> momentumHamiltonian(double J_ratio, int N);

void setHElement_momentum(double J_ratio, int N, std::list<Eigen::MatrixXcd> &H_subSubspace_list, int k,
                          const std::vector<int> &s_vector_k, const std::vector<int> &R_vector, int K);

std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXcd> & H_list, int N);

std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXd> & H_list, int N);

std::vector<double> getEnergiesFromBlocks(const std::list<std::list<Eigen::MatrixXcd>> & H_list, int N);

Eigen::MatrixXcd blkdiag(const std::list<Eigen::MatrixXcd>& matrix_list, int totalSize);

std::list<Eigen::MatrixXcd> blkdiag(const std::list<std::list<Eigen::MatrixXcd>> & matrix_doubleList);

void saveEnergies(const std::vector<double> & ergs, const std::string & path);


#endif //SPINCHAINED_HAMILTONIANBUILDERS_H
