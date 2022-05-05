//
// Created by mmaschke on 05/05/22.
//

#ifndef SPINCHAINED_DATAGENERATORS_H
#define SPINCHAINED_DATAGENERATORS_H

#include <utility>
#include <thread>
#include <mutex>
#include <algorithm>
#include <cmath>
#include <list>
#include <vector>
#include <complex>

#include <Eigen/Dense>

#include "hamiltonianBuilders.h"
#include "thermodynamics.h"

void saveExcitationErgsForVaryingJ(int N, int dataPointNum, double start, double end, std::string path);

void saveSpecificHeatForVaryingJ(int N, int dataPointNum, double betaOrT, double start, double end,
                                 bool isBeta, std::string path);

void saveSpecificHeatsForVaryingTemp(int N, int dataPointNum, double J_ratio, double start, double end,
                                     bool isBeta, std::string path);

void saveSusceptibilitiesForVaryingJ(int N, int dataPointNum, double betaOrT, double start, double end,
                                     bool isBeta, std::string path);

void saveSusceptibilitesForVaryingTemp(int N, int dataPointNum, double J_ratio, double start, double end,
                                       bool isBeta, std::string path);

template <typename T, typename U>
void savePairsToFile(std::list<std::pair<T, U>> pairList, std::string path);

std::vector<double> getMomentumErgsThreaded(const std::list<std::list<Eigen::MatrixXcd>> & H_list, int N);

std::vector<std::vector<double>> diagonalizeThreaded(const std::vector<double> & J_ratios, int N);

void writeThreadSafe (std::vector<std::vector<double>> & writeTo, const std::vector<double> & writeFrom);

void writeThreadSafe (std::vector<double> & writeTo, const std::vector<double> & writeFrom);



#endif //SPINCHAINED_DATAGENERATORS_H
