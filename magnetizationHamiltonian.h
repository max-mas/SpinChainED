//
// Created by MaxM on 13/05/2022.
//

#ifndef SPINCHAINED_MAGNETIZATIONHAMILTONIAN_H
#define SPINCHAINED_MAGNETIZATIONHAMILTONIAN_H

#include <vector>
#include <list>

#include <Eigen/Dense>

#include "helperFunctions.h"

// For a given m, get the corresponding block of the Hamiltonian.
Eigen::MatrixXd getMagnetizationBlock(double J_ratio, double m, int N);

// Get a list of all magnetization-blocks for a given N.
std::list<Eigen::MatrixXd> magnetizationHamiltonian(double J_ratio, int N);

// Set matrix element for a given state (magnetization approach).
void setHElement_magnetization(double J_ratio, int N, Eigen::MatrixXd & H_subspace_block,
                               const std::vector<int> &s_vector_m, int k);

#endif //SPINCHAINED_MAGNETIZATIONHAMILTONIAN_H
