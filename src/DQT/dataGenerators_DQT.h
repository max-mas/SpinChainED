#ifndef SPINCHAINED_DATAGENERATORS_DQT_H
#define SPINCHAINED_DATAGENERATORS_DQT_H

#include <string>
#include <random>
#include <utility>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "momentumHamiltonian_sparse.h"
#include "naiveHamiltonian_sparse.h"
#include "matrixExpIteration.h"

#include "../ED/dataGenerators.h"

void saveSpecificHeatsForVaryingTemp_DQT_parallel(int N, int dataPointNum, double J_ratio, double end, std::string path);

void saveSpecificHeatsForVaryingTemp_DQT_avg(int N, int dataPointNum, double J_ratio, double end, const std::string & path, int numOfRuns);

void saveSusceptibilityForVaryingTemp_DQT_avg(int N, int dataPointNum, double J_ratio, double end, const std::string & path, int numOfRuns);

void calcAbsDQTError(const std::string & EDpath, const std::string & DQTpath, const std::string & outPath);

std::vector<std::pair<double, double>> readPairVectorFromFile(const std::string & path);

void normaliseListOfVectors(std::vector<Eigen::VectorXcd> & vec);

Eigen::VectorXcd randomComplexVectorNormalised(int vecSize, double stdDev);

#endif //SPINCHAINED_DATAGENERATORS_DQT_H
