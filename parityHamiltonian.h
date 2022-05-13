//
// Created by MaxM on 13/05/2022.
//

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




std::list<std::list<std::list<Eigen::MatrixXd>>> parityHamiltonian(double J_ratio, int N);

double h_Element_parity(int a, int b, double l, double q, double k, double p, double N,
                        const std::vector<int> & s_vec, const std::vector<int> & R_vec, const std::vector<int> & m_vec);

double g_k(double k, double N);

double N_a_sigma(double g, double N, double sigmaR, double p, double k, double m);

double E_z_parity(int s, double J_ratio, int N);

#endif //SPINCHAINED_PARITYHAMILTONIAN_H
