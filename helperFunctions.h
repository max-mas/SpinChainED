//
// Created by mmaschke on 26/04/22.
//

#ifndef HEISENBERG_CHAIN_1D_C__HELPERFUNCTIONS_H
#define HEISENBERG_CHAIN_1D_C__HELPERFUNCTIONS_H

#include <vector>
#include <algorithm>

#include <cmath>
#include <complex>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "constants.h"

bool getBit(const int s, const int i);

void flipBit(int & s, int i);

int bitSum(const int s);

int findState(const std::vector<int> &s_list, const int s);

void cycleBits(int & s, const int N);

void setBit(int & s, const int i, const bool val);

int checkState(const int s, const int k, const int N);

void reflectBits(int & s, int N);

std::vector<int> checkState_parity(int s, int k, int N);

std::vector<int> representative(const int s, const int N);

std::vector<int> representative_parity(const int s, const int N);

double E_z(int s, int N);

double N_sigma(int s, int N, int R, int k, int sigma, int m, int p);

double hElement(int a, int b, int l, int q, int k, int p, int sigma);

void printMatrix(const Eigen::MatrixXd & M);

void printEnergies(const Eigen::VectorXd & v);
void printEnergies(const std::vector<double> & v);

#endif //HEISENBERG_CHAIN_1D_C__HELPERFUNCTIONS_H
