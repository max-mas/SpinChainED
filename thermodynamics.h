//
// Created by mmaschke on 03/05/22.
//

#ifndef SPINCHAINED_THERMODYNAMICS_H
#define SPINCHAINED_THERMODYNAMICS_H

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <complex>

double partitionFunction(const std::vector<double> & ergs, double beta);

double specificHeat(const std::vector<double> & ergs, double beta);

double susceptibility(const std::vector<double> & ergs, double beta, const Eigen::MatrixXcd & U,
                      const Eigen::MatrixXd & S_2);


#endif //SPINCHAINED_THERMODYNAMICS_H
