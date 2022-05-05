//
// Created by mmaschke on 03/05/22.
//

#ifndef SPINCHAINED_THERMODYNAMICS_H
#define SPINCHAINED_THERMODYNAMICS_H

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <complex>

double partitionFunction(const std::vector<double> & ergs, double beta, bool isBeta);

double specificHeat(const std::vector<double> & ergs, double betaOrT, bool isBeta) ;

double susceptibilityOld(const std::vector<double> & ergs, double betaOrT, bool isBeta, const Eigen::MatrixXcd & U,const Eigen::MatrixXd & S_2);

double susceptibility(const std::vector<double> & ergs, double betaOrT, bool isBeta, const Eigen::MatrixXcd & U,const Eigen::MatrixXd & S_2);


#endif //SPINCHAINED_THERMODYNAMICS_H
