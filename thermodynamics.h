//
// Created by mmaschke on 03/05/22.
//

#ifndef SPINCHAINED_THERMODYNAMICS_H
#define SPINCHAINED_THERMODYNAMICS_H

#include <vector>
#include <cmath>

double partitionFunction(const std::vector<double> & ergs, double beta);

double specificHeat(const std::vector<double> & ergs, double beta);

#endif //SPINCHAINED_THERMODYNAMICS_H
