//
// Created by mmaschke on 03/05/22.
//

#include "thermodynamics.h"

using std::vector;

double partitionFunction(const vector<double> & ergs, double beta) {
    double Z = 0;
    for (double erg : ergs) {
        Z += std::exp( -beta * erg );
    }
    return Z;
}

double specificHeat(const vector<double> & ergs, double beta) {
    double Z = partitionFunction(ergs, beta);
    double avg_H_sq = 0;
    double avg_H = 0;
    for (double erg : ergs) {
        avg_H_sq += std::exp(-beta * erg) * pow(erg, 2);
        avg_H    += std::exp(-beta * erg) * erg;
    }
    avg_H_sq *= 1/Z;
    avg_H    *= 1/Z;
    return pow(beta,2) * ( avg_H_sq - pow(avg_H, 2) );
}


