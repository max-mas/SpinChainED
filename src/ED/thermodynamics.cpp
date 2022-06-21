//
// This cpp file contains functions used to calculate thermodynamic quantities (partition function, specific heat
// and susceptibility).
//

#include "thermodynamics.h"

using std::vector;

// Calculates the partition function at a given temperature/beta for a given vector of eigenenergies.
double partitionFunction(const vector<double> & ergs, double betaOrT, bool isBeta) {
    double Z = 0;
    double temp_representative;
    if (isBeta) {
        temp_representative = betaOrT;
    } else {
        temp_representative = 1/betaOrT;
    }
    for (double erg : ergs) {
        Z += std::exp( -temp_representative * erg );
    }
    return Z;
}

// Calculates the specific heat for a given vector of eigenenergies at a given temperature/beta.
// Note: This could be made more efficient by using only the m = 0 block (TODO)
double specificHeat(const vector<double> & ergs, double betaOrT, bool isBeta) {
    double Z = partitionFunction(ergs, betaOrT, isBeta);

    double temp_representative;
    if (isBeta) {
        temp_representative = betaOrT;
    } else {
        temp_representative = 1/betaOrT;
    }

    double avg_H_sq = 0;
    double avg_H = 0;
    for (double erg : ergs) {
        avg_H_sq += std::exp(-temp_representative * erg) * pow(erg, 2);
        avg_H    += std::exp(-temp_representative * erg) * erg;
    }
    avg_H_sq *= 1/Z;
    avg_H    *= 1/Z;
    return pow(temp_representative,2) * ( avg_H_sq - pow(avg_H, 2) );
}

// Calculates the magnetic susceptibility at a given temperature/beta using <m_z²> = <S²>. U is the transformation
// matrix containing all eigenvectors of H.
double susceptibility(const vector<double> & ergs, double betaOrT, bool isBeta, const Eigen::MatrixXcd & T) {
    double S_2_avg = 0;
    double Z = 0;

    double temp_representative;
    if (isBeta) {
        temp_representative = betaOrT;
    } else {
        temp_representative = 1/betaOrT;
    }

    for (int i = 0; i < T.cols(); i++) {
        double S = -0.5 + sqrt(0.25 + T(i, i).real());
        Z += std::exp( -temp_representative * ergs[i] ) * (2 * S + 1);
        S_2_avg += std::exp(-temp_representative * ergs[i]) * S * (S + 1) * (2 * S + 1);
    }

    S_2_avg = S_2_avg/(3.0*Z);
    return temp_representative * S_2_avg ;
}


