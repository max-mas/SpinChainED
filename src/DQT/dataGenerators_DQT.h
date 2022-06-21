/**
 *  This file contains methods used to calculate thermodynamic quantities using dynamic quantum typicality (DQT).
 */


#ifndef SPINCHAINED_DATAGENERATORS_DQT_H
#define SPINCHAINED_DATAGENERATORS_DQT_H

#include <string>
#include <random>
#include <utility>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unsupported/Eigen/MatrixFunctions>

#include "momentumHamiltonian_sparse.h"
#include "naiveHamiltonian_sparse.h"
#include "matrixExpIteration.h"

#include "../ED/dataGenerators.h"

// Block-parallelized approach for calculation of specific heat using DQT.
void saveSpecificHeatsForVaryingTemp_DQT_parallel(int N, int dataPointNum, double J_ratio, double end, std::string path);

// Calculates specific heat using DQT and average over a number of runs. Also calculates standard deviation.
void saveSpecificHeatsForVaryingTemp_DQT_avg(int N, int dataPointNum, double J_ratio, double end, const std::string & path, int numOfRuns);

// Calculates susceptibility heat using DQT and average over a number of runs.
void saveSusceptibilityForVaryingTemp_DQT_avg(const int N, const int dataPointNum, const double J_ratio,
                                              const double end,const Eigen::SparseMatrix<std::complex<double>> & S2,
                                              const std::string & path, const int numOfRuns);

// Calculates partition function using DQT and average over a number of runs. Also calculates standard deviation.
void savePartitionFunction_DQT(const int N, const int dataPointNum, const double J_ratio, const double end, const std::string & path, const int numOfRuns);

// Calculates absolute error using QT and ED data from files with equal number of entries and equal dBeta.
void calcAbsDQTError(const std::string & EDpath, const std::string & DQTpath, const std::string & outPath);

// Calculates relative error using QT and ED data from files with equal number of entries and equal dBeta.
void calcRelDQTError(const std::string & EDpath, const std::string & DQTpath, const std::string & outPath);

// Reads double pairs from a file.
std::vector<std::pair<double, double>> readPairVectorFromFile(const std::string & path);

// Normalises a list of complex vectors to be interpreted as one large vector.
void normaliseListOfVectors(std::vector<Eigen::VectorXcd> & vec);

// Generates a random complex vector using a true random seed for a pseudo random number generator.
Eigen::VectorXcd randomComplexVectorNormalised(int vecSize, double stdDev);

#endif //SPINCHAINED_DATAGENERATORS_DQT_H
