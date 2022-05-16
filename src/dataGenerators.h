/**
 * This cpp file contains functions used for generating data of interest (excitation energies, specific heat and
 * susceptibilities) for later plotting in Python.
*/

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

#include "parityHamiltonian.h"
#include "naiveHamiltonian.h"
#include "magnetizationHamiltonian.h"
#include "momentumHamiltonian.h"
#include "diagonalizationMethods.h"
#include "thermodynamics.h"

// Saves excitation energy (i.e. energy difference between ground and 1st excited state) vor varying values of J1/J2.
// Note: If ground state is degenerate, zero is returned.
void saveExcitationErgsForVaryingJ(int N, int dataPointNum, double start, double end, const std::string & path);

// Saves specific heat at given temperature/beta for varying values of J1/J2.
void saveSpecificHeatForVaryingJ(int N, int dataPointNum, double betaOrT, double start, double end,
                                 bool isBeta, std::string path);

// Saves specific heat for a given value of J1/J2 for varying temperature/beta.
void saveSpecificHeatsForVaryingTemp(int N, int dataPointNum, double J_ratio, double start, double end,
                                     bool isBeta, std::string path);

// Saves susceptibilities at a given temperature/beta for varying values of J1/J2.
void saveSusceptibilitiesForVaryingJ(int N, int dataPointNum, double betaOrT, double start, double end,
                                     bool isBeta, std::string path);

// Saves susceptibilities for a given value of J1/J2 for varying temperature/beta.
void saveSusceptibilitesForVaryingTemp(int N, int dataPointNum, double J_ratio, double start, double end,
                                       bool isBeta, std::string path);

// Saves tuples of (E, k) to file.
void saveEnergyDispersion(int N, double J_ratio, std::string path);

// Saves tuples of (E, m, k) to file.
void saveEnergyDispersionWithMag(int N, double J_ratio, std::string path);

// Used to write data tuples to a file at path. Operator << must be defined correctly.
template <typename T, typename U>
void savePairsToFile(const std::list<std::pair<T, U>> & pairList, std::string path);

template <typename T, typename U, typename V>
void saveTripleToFile(const std::list<std::tuple<T, U, V>> & pairList, std::string path);

// Enables threaded writing to a vector.
void writeThreadSafe (std::vector<std::vector<double>> & writeTo, const std::vector<double> & writeFrom);

// Enables threaded writing to a vector.
void writeThreadSafe (std::vector<double> & writeTo, const std::vector<double> & writeFrom);



#endif //SPINCHAINED_DATAGENERATORS_H