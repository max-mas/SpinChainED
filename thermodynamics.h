//
// This cpp file contains functions used to calculate thermodynamic quantities (partition function, specific heat
// and susceptibility).
//

#ifndef SPINCHAINED_THERMODYNAMICS_H
#define SPINCHAINED_THERMODYNAMICS_H

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <complex>

// Calculates the partition function at a given temperature/beta for a given vector of eigenenergies.
double partitionFunction(const std::vector<double> & ergs, double beta, bool isBeta);

// Calculates the specific heat for a given vector of eigenenergies at a given temperature/beta.
// Note: This could be made more efficient by using only the m = 0 block
double specificHeat(const std::vector<double> & ergs, double betaOrT, bool isBeta) ;

// Calculates the magnetic susceptibility at a given temperature/beta using <m_z²> = <S²>. U is the transformation
// matrix containing all eigenvectors of H.
double susceptibility(const std::vector<double> & ergs, double betaOrT, bool isBeta, const Eigen::MatrixXcd & U,const Eigen::MatrixXd & S_2);


#endif //SPINCHAINED_THERMODYNAMICS_H
