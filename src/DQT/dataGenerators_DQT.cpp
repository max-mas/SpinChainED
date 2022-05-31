#include "dataGenerators_DQT.h"

using std::string;
using std::complex;
using std::vector;
using std::list;

using Eigen::VectorXcd;
using Eigen::SparseMatrix;


void saveSpecificHeatsForVaryingTemp_DQT_parallel(int N, int dataPointNum, double J_ratio, double end, string path) {
    double beta = 0;
    double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<double> Cs;

     if (N % 2 == 0 && N >= 6) {
        const vector<SparseMatrix<complex<double>>> H  = momentumHamiltonian_sparse_blocks(J_ratio, N);
        vector<VectorXcd> psi_vec;

        for (const auto & h : H) {
            psi_vec.emplace_back(randomComplexVectorNormalised(h.cols(), 1));
        }
        normaliseListOfVectors(psi_vec);

        for (int i = 0; i < dataPointNum; i++) {
            vector<double> avg_H2_vec;
            vector<double> avg_H_vec;
#pragma omp parallel for default(none) shared(psi_vec, avg_H2_vec, avg_H_vec, dBeta, i)
            for (int j = 0; j < H.size(); j++) {
                double avg_H2_block = psi_vec[j].dot(H[j] * (H[j] * psi_vec[j])).real();
                double avg_H_block  = psi_vec[j].dot(H[j] * psi_vec[j]).real();
                writeThreadSafe(avg_H2_vec, {avg_H2_block});
                writeThreadSafe(avg_H_vec, {avg_H_block});
                iterateState_beta(H[j], psi_vec[j], dBeta);
            }
            double avg_H2 = std::accumulate(avg_H2_vec.begin(), avg_H2_vec.end(), 0.0);
            double avg_H  = std::accumulate(avg_H_vec.begin() , avg_H_vec.end() , 0.0);

            double diff = avg_H2 - pow(avg_H, 2);
            double C = pow(beta, 2) * ( diff ) / N;
            Cs.emplace_back(C);

            normaliseListOfVectors(psi_vec);

            beta += dBeta;
        }
    }

    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(betas[j], Cs[j]) );
    }
    savePairsToFile(out, std::move(path));
}

void saveSpecificHeatsForVaryingTemp_DQT_avg(const int N, const int dataPointNum, const double J_ratio, const double end, const string & path, const int numOfRuns) {
    const double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<double> Cs(dataPointNum, 0);

    if (N % 2 == 0 && N >= 6) {
        const SparseMatrix<complex<double>> H = momentumHamiltonian_sparse(J_ratio, N);

#pragma omp parallel for default(none) shared(Cs)
        for (int k = 0; k < numOfRuns; k++) {
            VectorXcd psi = randomComplexVectorNormalised(pow(2, N), 1);
            double beta = 0;
            for (int i = 0; i < dataPointNum; i++) {
                double avg_H2 = psi.dot(H * (H * psi)).real();
                double avg_H = psi.dot(H * psi).real();

                double diff = avg_H2 - pow(avg_H, 2);
                double C = pow(beta, 2) * (diff) / N;
#pragma omp critical
                Cs[i] += C;

                beta += dBeta;
                iterateState_beta(H, psi, dBeta);
                psi.normalize();
            }
        }

        for (double & C : Cs) {
            C /= (double) numOfRuns;
        }
    }

    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(betas[j], Cs[j]) );
    }
    savePairsToFile(out, path);
}

void saveSusceptibilityForVaryingTemp_DQT(const int N, const int dataPointNum, const double J_ratio, const double end, const string & path) {
    const double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<double> Xs(dataPointNum, 0);

    if (N % 2 == 0 && N >= 6) {
        const SparseMatrix<complex<double>> S2 = spinOp2_momentum_sparse(N);
        const SparseMatrix<complex<double>> H  = momentumHamiltonian_sparse(J_ratio, N);

        VectorXcd psi = randomComplexVectorNormalised(pow(2, N), 1);
        double beta = 0;
        for (int i = 0; i < dataPointNum; i++) {
            double avg_S2 = psi.dot(S2* psi).real();

            double x = beta * avg_S2 / (3.0 * N);

            Xs[i] += x;

            beta += dBeta;
            iterateState_beta(H, psi, dBeta);
            psi.normalize();
        }
    }

    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(betas[j], Xs[j]) );
    }
    savePairsToFile(out, path);
}

void normaliseListOfVectors(vector<VectorXcd> & vec) {
    double norm2 = 0;
    for (VectorXcd & v : vec) {
        norm2 += pow (v.norm(),2);
    }
    double norm = sqrt(norm2);
    for (VectorXcd & v : vec) {
        v /= norm;
    }
}

VectorXcd randomComplexVectorNormalised(int vecSize, double stdDev) {
    VectorXcd psi(vecSize);

    std::random_device generator; //may not be available
    std::mt19937 gen(generator());
    std::normal_distribution<double> distribution(0.0, stdDev);

    for (int i = 0; i < vecSize; i++) {
        double re = distribution(gen);
        gen();
        double im = distribution(gen);
        gen();
        psi(i) = complex<double>(re, im);
    }

    psi.normalize();

    return psi;
}
