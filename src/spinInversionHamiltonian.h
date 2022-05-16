//
// Created by mmaschke on 16/05/22.
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

std::list<std::list<std::list<Eigen::MatrixXd>>> spinInversionHamiltonian(double J_ratio, int N);

void setHElementsForState_inversion(int N, int k, int p, int z, const std::vector<int> & s_vector_k,
                                    const std::vector<int> & R_vector, const std::vector<int> &m_vector,
                                    const std::vector<int> & n_vector, const std::vector<int> & c_vector,
                                    int K, Eigen::MatrixXd & H, int a, const int s, int n, int neighbour, double J);

double h_Element_inversion(int a, int b, double l, double q, double g, double k, double p, double z, double N,
                           const std::vector<int> & s_vec, const std::vector<int> & R_vec, const std::vector<int> & m_vec,
                           const std::vector<int> & n_vec, const std::vector<int> & c_vec);

double N_a_sigma_inversion(double N, int c, double sigma, double R, double p, double k, double z, double m, double n);

void getStates_k_p_z(int N, const std::vector<int> &s_vector_m, int k, int p, int z, std::vector<int> &s_vector_k,
                     std::vector<int> &R_vector, std::vector<int> &m_vector, std::vector<int> & n_vector,
                     std::vector<int> & c_vector);

std::vector<int> checkState_inversion(int s, int k, int N);

std::vector<int> representative_inversion(const int s, const int N);

int getStateClass(int m_p, int m_z, int m_zp);

void invertBits(int & s, int N);

#endif //SPINCHAINED_SPININVERSIONHAMILTONIAN_H
