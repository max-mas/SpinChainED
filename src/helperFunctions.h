/**
    This cpp file contains functions used to manipulate integers in the way needed for the chosen numerical
    representation of spin states, as well as functions used for console output.
*/

#ifndef HEISENBERG_CHAIN_1D_C__HELPERFUNCTIONS_H
#define HEISENBERG_CHAIN_1D_C__HELPERFUNCTIONS_H

#include <vector>
#include <algorithm>

#include <cmath>
#include <complex>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "constants.h"

// Gets bit in integer s at index i.
bool getBit(const int s, const int i);

// Flips bit in int s at index i using bitwise XOR.
void flipBit(int & s, int i);

// Sets bit in int s at index i using bitwise AND.
void setBit(int & s, const int i, const bool val);

// Counts set (=1) bits in int s using bitwise AND.
int bitSum(const int s);

// Returns index of first occurrence of int s in s_list. Returns -1 if not found.
int findState(const std::vector<int> &s_list, const int s);

// Cyclic left bit shift by 1 for integers in which only the first N bits must be set! Otherwise, undefined behavior.
void cycleBits(int & s, const int N);

// Calls cycleBits twice.
void cycleBits2(int &s, const int N);

// Get vector containing all states corresponding to a given m.
std::vector<int> getStates_m(int N, int n_up);

// Reflects bits about the center of the chain.
void reflectBits(int & s, int N);

// Prints real-valued matrix to console.
void printMatrix(const Eigen::MatrixXd & M);
// Prints complex-valued matrix to console.
void printMatrix(const Eigen::MatrixXcd & M);

// Prints doubles to console from Eigen::Vector.
void printEnergies(const Eigen::VectorXd & v);
// Prints doubles to console from std::vector.
void printEnergies(const std::vector<double> & v);

long fact(int n);

#endif //HEISENBERG_CHAIN_1D_C__HELPERFUNCTIONS_H
