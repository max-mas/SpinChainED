//
// Created by mmaschke on 4/25/22.
//

#ifndef HEISENBERG_CHAIN_1D_C_MAIN_H
#define HEISENBERG_CHAIN_1D_C_MAIN_H

#include "helperFunctions.h"
#include "hamiltonianBuilders.h"
#include "thermodynamics.h"

std::vector<double> getMomentumErgsThreaded(const std::list<std::list<Eigen::MatrixXcd>> & H_list, int N);

std::vector<std::vector<double>> diagonalizeThreaded(const std::vector<double> & J_ratios, int N);

void writeThreadSafe (std::vector<std::vector<double>> & writeTo, const std::vector<double> & writeFrom);

void writeThreadSafe (std::vector<double> & writeTo, const std::vector<double> & writeFrom);

#endif //HEISENBERG_CHAIN_1D_C_MAIN_H
