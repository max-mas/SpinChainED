/**
 * This file contains methods used to construct the (m, k, p) H-blocks using magnetization conservation as well as
 * translational and mirror lattice symmetries.
 */

#ifndef SPINCHAINED_PARITYHAMILTONIAN_H
#define SPINCHAINED_PARITYHAMILTONIAN_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <list>

#include <stdexcept>
#include <string>

#include <cmath>
#include <complex>

#include <Eigen/Dense>

#include "helperFunctions.h"

// Returns (m, k, p) blocks for given N and J1/J2.
std::list<std::list<std::list<Eigen::MatrixXd>>> parityHamiltonian(double J_ratio, int N);

// Returns off-diagonal operator matrix elements.
double h_Element_parity(int a, int b, double l, double q, double k, double p, double N,
                        const std::vector<int> & s_vec, const std::vector<int> & R_vec, const std::vector<int> & m_vec);

double g_k(double k, double N);

// Returns state normalization constant.
double N_a_sigma(double g, double N, double sigmaR, double p, double k, double m);

// Returns diagonal operator matrix elements.
double E_z_parity(int s, double J_ratio, int N);

// Checks whether state is valid representative and returns flags if not so, else periodicities
// (with reflection and without) are returned.
std::vector<int> checkState_parity(const int s, const int k, const int N);

// Determines the correct representative of state s.
std::vector<int> representative_parity(const int s, const int N);

#endif //SPINCHAINED_PARITYHAMILTONIAN_H
