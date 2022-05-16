//
// Created by MaxM on 16/05/2022.
//

#ifndef SPINCHAINED_SPININVERSION2_H
#define SPINCHAINED_SPININVERSION2_H

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

std::list<std::list<std::list<Eigen::MatrixXd>>> spinInversionHamiltonian2(double J_ratio, int N);

void setHElementsForState_inversion2(int N, int k, int p, int z, const std::vector<int> & s_vec,
                                    const std::vector<int> & R_vec, const std::vector<int> & m_vec,
                                    const std::vector<int> & n_vec, const std::vector<int> & c_vec, int K,
                                    Eigen::MatrixXd & H, int a, int s, int n, int neighbour, double J);

double hElement_inversion2(const int a, const int b, const double l, const double q, const double g, const double k,
                          const double p, const double z, const double N, const std::vector<int> & R_vec,
                          const std::vector<int> & m_vec, const std::vector<int> & n_vec,
                          const std::vector<int> & c_vec);

double N_a_sigma_inv2(const double N, const double sigmaR, const double k, const double p, const double z,
                      const double m, const double n, const int c);

int getStateCase2(const int m_p, const int m_z, const int m_pz);

std::vector<int> checkState_inversion2(const int s, const int k, const int N);

std::vector<int> representative_inversion2(const int s, const int N);

void invertBits2(int & s, int N);

#endif //SPINCHAINED_SPININVERSION2_H
