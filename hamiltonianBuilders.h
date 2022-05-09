//
// This cpp file contains functions used to generate the Hamiltonian using different approaches, as well as the spin
// operator S². Also included: Functions to get eigenvalues and generate full matrices from blocks.
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

// Generate full-size S² for given N.
Eigen::MatrixXd spinOperator_sq(int N);

// Generate S² for a given subset of states (in practice: m = 0 states).
Eigen::MatrixXd spinOperator_sq(std::vector<int> states, int N);

// Generate full-size hamiltonian using the naive approach.
Eigen::MatrixXd naiveHamiltonian(double J_ratio, int N);

// Set matrix element for a given state (naive approach).
void setHElement_naive(double J_ratio, int N, Eigen::MatrixXd &H, int a);

// For a given m, get the corresponding block of the Hamiltonian.
Eigen::MatrixXd getMagnetizationBlock(double J_ratio, double m, int N);

// Get a list of all magnetization-blocks for a given N.
std::list<Eigen::MatrixXd> magnetizationHamiltonian(double J_ratio, int N);

// Get vector containing all states corresponding to a given m.
std::vector<int> getStates_m(int N, int n_up);

// Set matrix element for a given state (magnetization approach).
void setHElement_magnetization(double J_ratio, int N, Eigen::MatrixXd & H_subspace_block,
                               const std::vector<int> &s_vector_m, int k);

// Generate list of list of (m, k)-blocks of the Hamiltonian for a given N.
std::list<std::list<Eigen::MatrixXcd>> momentumHamiltonian(double J_ratio, int N);

// Set matrix element for a given state (momentum-approach).
void setHElement_momentum(double J_ratio, int N, std::list<Eigen::MatrixXcd> &H_subSubspace_list, int k,
                          const std::vector<int> &s_vector_k, const std::vector<int> &R_vector, int K);

// Diagonalization-methods for a number of cases.
std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXcd> & H_list, int N);
std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXcd> & H_list);
std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXd> & H_list);
std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXd> & H_list, int N);
std::vector<double> getEnergiesFromBlocks(const std::list<std::list<Eigen::MatrixXcd>> & H_list, int N);

std::vector<std::vector<std::vector<double>>> getEnergiesFromBlocksByK(
        const std::list<std::list<Eigen::MatrixXcd>> & H_list);

// Generate full-sized matrix from blocks.
Eigen::MatrixXcd blkdiag(const std::list<Eigen::MatrixXcd>& matrix_list, int totalSize);
std::list<Eigen::MatrixXcd> blkdiag(const std::list<std::list<Eigen::MatrixXcd>> & matrix_doubleList);

// Save vector containing eigenvalues to file.
void saveEnergies(const std::vector<double> & ergs, const std::string & path);


#endif //SPINCHAINED_HAMILTONIANBUILDERS_H
