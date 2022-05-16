/**
 * This cpp file contains methods used to get Eigenvalues from lists of self-adjoint H-blocks.
*/

#ifndef SPINCHAINED_DIAGONALIZATIONMETHODS_H
#define SPINCHAINED_DIAGONALIZATIONMETHODS_H

#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include <complex>

#include <Eigen/Dense>

#include "dataGenerators.h"

// Diagonalization-methods for a number of cases.
std::vector<double> getEnergiesFromBlocks(const std::list<std::list<std::list<Eigen::MatrixXd>>> & H_list);
std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXcd> & H_list);
std::vector<double> getEnergiesFromBlocks(const std::list<Eigen::MatrixXd> & H_list);
std::vector<double> getEnergiesFromBlocks(const std::list<std::list<Eigen::MatrixXcd>> & H_list);
std::vector<double> getEnergiesFromBlocks(const std::list<std::list<Eigen::MatrixXd>> & H_list);

// Returns energies sorted by m and k, used for generation of dispersion plots.
std::vector<std::vector<std::vector<double>>> getEnergiesFromBlocksByK(
        const std::list<std::list<Eigen::MatrixXcd>> & H_list);

// Threaded version of getEnergiesFromBlocks for the momentum state ansatz.
std::vector<double> getMomentumErgsThreaded(const std::list<std::list<Eigen::MatrixXcd>> & H_list, int N);

// Threaded version of getEnergiesFromBlocks for the semi-momentum/parity state ansatz.
std::vector<double> getParityErgsThreaded(const std::list<std::list<std::list<Eigen::MatrixXd>>> & H_list, int N);

// Returns all eigenvalues for a given vector of values of J1/J2. Threaded.
std::vector<std::vector<double>> diagonalizeThreaded(const std::vector<double> & J_ratios, int N);

// Generate full-sized matrix from blocks.
Eigen::MatrixXcd blkdiag(const std::list<Eigen::MatrixXcd>& matrix_list, int totalSize);
std::list<Eigen::MatrixXcd> blkdiag(const std::list<std::list<Eigen::MatrixXcd>> & matrix_doubleList);

// Save vector containing eigenvalues to file.
void saveEnergies(const std::vector<double> & ergs, const std::string & path);

#endif //SPINCHAINED_DIAGONALIZATIONMETHODS_H
