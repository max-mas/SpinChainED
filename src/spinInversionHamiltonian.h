//
// Created by MaxM on 16/05/2022.
//

#ifndef SPINCHAINED_SPININVERSIONHAMILTONIAN_H
#define SPINCHAINED_SPININVERSIONHAMILTONIAN_H

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
#include "parityHamiltonian.h"

std::list<std::list<std::list<Eigen::MatrixXd>>> spinInversionHamiltonian(double J_ratio, int N, int n_up_min, int n_up_max);

std::vector<double> getEnergies_memorySaving_threaded_inversion(double J_ratio, int N);

void getStates_k_p_z(int N, std::vector<int> &s_vector_m, int k, int p, int z, std::vector<int> &s_vector_k,
                     std::vector<int> &R_vector, std::vector<int> &m_vector, std::vector<int> &n_vector,
                     std::vector<int> &c_vector);

void setHElementsForState_inversion(int N, int k, int p, int z, const std::vector<int> & s_vec,
                                    const std::vector<int> & R_vec, const std::vector<int> & m_vec,
                                    const std::vector<int> & n_vec, const std::vector<int> & c_vec, int K,
                                    Eigen::MatrixXd & H, int a, int s, int n, int neighbour, double J);

double hElement_inversion(const int a, const int b, const double l, const double q, const double g, const double k,
                          const double p, const double z, const double N, const std::vector<int> & R_vec,
                          const std::vector<int> & m_vec, const std::vector<int> & n_vec,
                          const std::vector<int> & c_vec);

double N_a_sigma_inv(const double N, const double sigmaR, const double k, const double p, const double z,
                     const double m, const double n, const int c);

int getStateClass(const int m_p, const int m_z, const int m_pz);

std::vector<int> checkState_inversion(const int s, const int k, const int N);

std::vector<int> representative_inversion(const int s, const int N);

void invertBits(int & s, int N);

#endif //SPINCHAINED_SPININVERSIONHAMILTONIAN_H
