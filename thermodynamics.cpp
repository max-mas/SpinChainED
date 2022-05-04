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

double susceptibility(const vector<double> & ergs, double beta, const Eigen::MatrixXcd & U,const Eigen::MatrixXd & S_2) {
    std::complex<double> S_2_avg = 0;
    Eigen::MatrixXcd S_2_transform =  U.inverse() * S_2 * U;
    double Z = partitionFunction(ergs, beta);

    for (double erg : ergs) {
        for (int i = 0; i < S_2.cols(); i++) {
            S_2_avg += std::exp(-beta * erg) * S_2_transform(i, i);
        }
    }
    double S_2_avg_real = (S_2_avg / 3.0).real();
    return beta / Z * S_2_avg_real;
}


