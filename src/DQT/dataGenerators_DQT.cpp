#include "dataGenerators_DQT.h"

using std::string;
using std::complex;
using std::vector;
using std::list;

using Eigen::VectorXcd;
using Eigen::SparseMatrix;


void saveSpecificHeatsForVaryingTemp_DQT(int N, int dataPointNum, double J_ratio, double end, string path) {
    double stdDev = 1;// exp(1);
    double beta = 0;
    double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<double> Cs;

     if (N % 2 == 0 && N >= 6) {
        SparseMatrix<complex<double>> H  = naiveHamiltonian_sparse(J_ratio, N).cast<complex<double>>();
        SparseMatrix<complex<double>> H2 = H*H;

        Eigen::VectorXcd psi = randomComplexVectorNormalised(N, stdDev);

        for (int i = 0; i < dataPointNum; i++) {
            complex<double> avg_H2 = (psi.adjoint() * H2 * psi)(0,0);
            complex<double> avg_H  = (psi.adjoint() * H  * psi)(0,0);

            double C = pow(beta, 2) * ( avg_H2.real() - pow(avg_H.real(), 2) ) / N;

            Cs.emplace_back(C);

            beta += dBeta;
            iterateState_beta(H, psi, dBeta);
        }
    }

    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(betas[j], Cs[j]) );
    }
    savePairsToFile(out, path);
}

VectorXcd randomComplexVectorNormalised(int N, double stdDev) {
    int vecSize = pow(2, N);
    VectorXcd psi(vecSize);

    std::random_device generator; //may not be available
    std::normal_distribution<double> distribution(0.0, stdDev);

    for (int i = 0; i < vecSize; i++) {
        double re = distribution(generator);
        double im = distribution(generator);
        psi(i) = complex<double>(re, im);
    }

    psi.normalize();

    return psi;
}
