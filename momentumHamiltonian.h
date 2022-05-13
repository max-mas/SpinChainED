/**
 * This file contains methods used to construct (m, k) H-blocks using momentum conservation and translational lattice
 * symmetry.
 */

#ifndef SPINCHAINED_MOMENTUMHAMILTONIAN_H
#define SPINCHAINED_MOMENTUMHAMILTONIAN_H

#include <vector>
#include <list>
#include <complex>

#include <Eigen/Dense>

#include "helperFunctions.h"

// Generate list of list of (m, k)-blocks of the Hamiltonian for a given N.
std::list<std::list<Eigen::MatrixXcd>> momentumHamiltonian(double J_ratio, int N);

// Set matrix element for a given state (momentum-approach).
void setHElement_momentum(double J_ratio, int N, std::list<Eigen::MatrixXcd> &H_subSubspace_list, int k,
                          const std::vector<int> &s_vector_k, const std::vector<int> &R_vector, int K);

// Checks if state int s is smallest among its translations and compatible with momentum k.
int checkState(const int s, const int k, const int N);

// Generates representative for given state int s. 1st return: representative. 2nd return:
// Needed number of translations.
std::vector<int> representative(const int s, const int N);



#endif //SPINCHAINED_MOMENTUMHAMILTONIAN_H
