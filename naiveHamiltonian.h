//
// Created by MaxM on 13/05/2022.
//

#ifndef SPINCHAINED_NAIVEHAMILTONIAN_H
#define SPINCHAINED_NAIVEHAMILTONIAN_H

#include <list>
#include <vector>

#include <Eigen/Dense>

#include "helperFunctions.h"

// Generate full-size S² for given N.
Eigen::MatrixXd spinOperator_sq(int N);

// Generate S² for a given subset of states (in practice: m = 0 states).
Eigen::MatrixXd spinOperator_sq(std::vector<int> states, int N);

// Generate full-size hamiltonian using the naive approach.
Eigen::MatrixXd naiveHamiltonian(double J_ratio, int N);

// Set matrix element for a given state (naive approach).
void setHElement_naive(double J_ratio, int N, Eigen::MatrixXd &H, int a);


#endif //SPINCHAINED_NAIVEHAMILTONIAN_H
