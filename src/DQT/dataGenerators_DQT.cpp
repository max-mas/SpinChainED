/**
 *  This file contains methods used to calculate thermodynamic quantities using dynamic quantum typicality (DQT).
 */

#include "dataGenerators_DQT.h"
#include <chrono>

using std::string;
using std::complex;
using std::vector;
using std::list;

using Eigen::VectorXcd;
using Eigen::SparseMatrix;

// Block-parallelized approach for calculation of specific heat using DQT.
void saveSpecificHeatsForVaryingTemp_DQT_parallel(int N, int dataPointNum, double J_ratio, double end, string path) {
    double beta = 0; // iteration always starts at beta = 0
    double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<double> Cs;

     if (N % 2 == 0 && N >= 6) {
        const vector<SparseMatrix<complex<double>>> H  = momentumHamiltonian_sparse_blocks(J_ratio, N, 0, N);
        vector<VectorXcd> psi_vec;

        for (const auto & h : H) {
            psi_vec.emplace_back(randomComplexVector(h.cols(), 1));
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

// Calculates specific heat using DQT and average over a number of runs. Also calculates standard deviation.
void
saveSpecificHeatsForVaryingTemp_DQT_avg(const int N, const int dataPointNum, const double J_ratio, const double end,
                                        const string &path, const int numOfRuns, bool bench) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();


    const double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<vector<double>> Cs(dataPointNum);
    vector<double> actualCs;
    vector<double> stdDevs;

    if (N % 2 == 0 && N >= 6) {
        const SparseMatrix<complex<double>> H = momentumHamiltonian_sparse(J_ratio, N);

#pragma omp parallel for default(none) shared(Cs)
        for (int k = 0; k < numOfRuns; k++) {
            VectorXcd psi = randomComplexVectorNormalised((int) pow(2, N), 1.0);
            double beta = 0;
            for (int i = 0; i < dataPointNum; i++) {
                double avg_H2 = psi.dot(H * (H * psi)).real();
                double avg_H = psi.dot(H * psi).real();

                double diff = avg_H2 - pow(avg_H, 2);
                double C = pow(beta, 2) * (diff) / N;
#pragma omp critical
                Cs[i].emplace_back(C);

                beta += dBeta;
                iterateState_beta(H, psi, dBeta);
                psi.normalize();
            }
        }

        for (vector<double> & C : Cs) {
            double mean = std::accumulate(C.begin(), C.end(), 0.0) / numOfRuns;
            double stdDev = 0.0;
            for (double c : C) {
                stdDev += pow( (mean - c), 2);
            }
            stdDev = sqrt(stdDev/(double) numOfRuns);
            actualCs.emplace_back(mean);
            stdDevs.emplace_back(stdDev);
        }
    }

    list<std::tuple<double, double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        std::tuple<double, double, double> a(betas[j], actualCs[j], stdDevs[j]);
        out.emplace_back( a );
    }
    saveTripleToFile(out, path);

    if (bench) {
        std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
        std::ifstream statFile("/proc/self/stat");
        std::string statLine;
        std::getline(statFile, statLine);
        std::istringstream iss(statLine);
        std::string entry;
        long long memUsage;
        for (int i = 1; i <= 24; i++) {
            std::getline(iss, entry, ' ');
            if (i == 24) {
                memUsage = stoi(entry);
            }
        }

        long time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count();
        std::cout << "N = " << N << ", time = " << time_ns/1e9 << " s, ram = " << 4096*memUsage/1e9 << " GB" << std::endl;
    }
}

// Calculates susceptibility heat using DQT and average over a number of runs.
void saveSusceptibilityForVaryingTemp_DQT_avg(const int N, const int dataPointNum, const double J_ratio,
                                              const double end, const SparseMatrix<complex<double>> & S2,
                                              const string & path, const int numOfRuns) {
    const double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<vector<double>> Xs(dataPointNum);
    vector<double> actualXs;
    vector<double> stdDevs;

    if (N % 2 == 0 && N >= 6) {
        const SparseMatrix<complex<double>> H  = momentumHamiltonian_sparse(J_ratio, N);

#pragma omp parallel for default(none) shared(Xs, S2, std::cout)
        for (int k = 0; k < numOfRuns; k++) {

            VectorXcd psi = randomComplexVectorNormalised((int) pow(2, N), 1.0);

            double beta = 0;
            for (int i = 0; i < dataPointNum; i++) {
                double avg_S2 = (psi.adjoint() * S2 * psi)(0, 0).real();

                double x = beta * avg_S2 / (3.0 * N);
                std::cout << avg_S2 << std::endl;

#pragma omp critical
                Xs[i].emplace_back(x);

                beta += dBeta;
                iterateState_beta(H, psi, dBeta);
                psi.normalize();
            }
        }

        for (vector<double> & X : Xs) {
            double mean = std::accumulate(X.begin(), X.end(), 0.0) / numOfRuns;
            double stdDev = 0.0;
            for (double x : X) {
                stdDev += pow( (mean - x), 2);
            }
            stdDev = sqrt(stdDev/(double) numOfRuns);
            actualXs.emplace_back(mean);
            stdDevs.emplace_back(stdDev);
        }
    }

    list<std::tuple<double, double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        std::tuple<double, double, double> a(betas[j], actualXs[j], stdDevs[j]);
        out.emplace_back( a );
    }
    saveTripleToFile(out, path);
}

// Calculates susceptibility heat using DQT and average over a number of runs.
void saveSusceptibilityForVaryingTemp_DQT_parallel(const int N, const int dataPointNum, const double J_ratio,
                                              const double end, const vector<SparseMatrix<double>> & S2_list,
                                              const string & path) {
    double beta = 0; // iteration always starts at beta = 0
    const double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<double> Xs;

    if (N % 2 == 0 && N >= 6) {
        const vector<SparseMatrix<double>> H  = magnetizationHamiltonian_sparse(J_ratio, N);
        vector<VectorXcd> psi_vec;

        for (const auto & h : H) {
            psi_vec.emplace_back(randomComplexVector(h.cols(), 1));
        }
        normaliseListOfVectors(psi_vec);

        for (int i = 0; i < dataPointNum; i++) {
            vector<double> avg_S2_vec;
#pragma omp parallel for default(none) shared(psi_vec, avg_S2_vec, i, S2_list)
            for (int j = 0; j < H.size(); j++) {
                double S2_avg_block = (psi_vec[j].adjoint() * (S2_list[j] * psi_vec[j]))(0, 0).real();

                writeThreadSafe(avg_S2_vec, {S2_avg_block});
                iterateState_beta(H[j], psi_vec[j], dBeta);
            }
            double avg_S2 = std::accumulate(avg_S2_vec.begin(), avg_S2_vec.end(), 0.0);
            double X = beta * avg_S2 /3.0/N;
            Xs.emplace_back(X);

            //std::cout << avg_S2 << std::endl;

            normaliseListOfVectors(psi_vec);


            beta += dBeta;
        }
    }

    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        std::pair<double, double> a(betas[j], Xs[j]);
        out.emplace_back( a );
    }
    savePairsToFile(out, path);
}

// Calculates partition function using DQT and average over a number of runs. Also calculates standard deviation.
void savePartitionFunction_DQT(const int N, const int dataPointNum, const double J_ratio, const double end,
                               const string & path, const int numOfRuns) {
    const double dBeta = end / (double) dataPointNum;
    Eigen::VectorXd betas = Eigen::VectorXd::LinSpaced(dataPointNum, 0, end);

    vector<vector<double>> Zs(dataPointNum);
    vector<double> actualZs;
    vector<double> stdDevs;

    if (N % 2 == 0 && N >= 6) {
        const Eigen::MatrixXcd H = momentumHamiltonian_sparse(J_ratio, N);
        const SparseMatrix<complex<double>> H_s = momentumHamiltonian_sparse(J_ratio, N);
        //const SparseMatrix<complex<double>> H2 = H*H;

#pragma omp parallel for default(none) shared(Zs)
        for (int k = 0; k < numOfRuns; k++) {
            VectorXcd psi = randomComplexVectorNormalised((int) pow(2, N), 1.0);
            double beta = 0;
            for (int i = 0; i < dataPointNum; i++) {
                //double avg_H2 = psi.dot(H * (H * psi)).real();
                //double avg_H = psi.dot(H * psi).real();

                //double Z = 1 - beta * avg_H + pow(beta, 2)/2.0 * avg_H2;
                double Z = (psi.dot( (-beta*H).exp() * psi)).real();
#pragma omp critical
                Zs[i].emplace_back(Z);

                beta += dBeta;
                iterateState_beta(H_s, psi, dBeta);
                psi.normalize();
            }
        }

        for (vector<double> & Z : Zs) {
            double mean = std::accumulate(Z.begin(), Z.end(), 0.0) / numOfRuns;
            double stdDev = 0.0;
            for (double z : Z) {
                stdDev += pow( (mean - z), 2);
            }
            stdDev = sqrt(stdDev/(double) numOfRuns);
            actualZs.emplace_back(mean);
            stdDevs.emplace_back(stdDev);
        }
    }

    list<std::tuple<double, double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        std::tuple<double, double, double> a(betas[j], actualZs[j], stdDevs[j]);
        out.emplace_back( a );
    }
    saveTripleToFile(out, path);
}

// Calculates absolute error using QT and ED data from files with equal number of entries and equal dBeta.
void calcAbsDQTError(const string & EDpath, const string & DQTpath, const string & outPath) {
    vector<std::pair<double, double>> EDData  = readPairVectorFromFile(EDpath );
    vector<std::pair<double, double>> DQTData = readPairVectorFromFile(DQTpath);
    list<std::pair<double, double>> betasAndDiffs;
    for (int i = 0; i < EDData.size(); i++) {
        betasAndDiffs.emplace_back(std::pair<double, double>( EDData[i].first, abs( EDData[i].second - DQTData[i].second ) ));
    }
    savePairsToFile(betasAndDiffs, outPath);
}

// Calculates relative error using QT and ED data from files with equal number of entries and equal dBeta.
// If value is too small (1e-10), it is ignored.
void calcRelDQTError(const string & EDpath, const string & DQTpath, const string & outPath) {
    vector<std::pair<double, double>> EDData  = readPairVectorFromFile(EDpath );
    vector<std::pair<double, double>> DQTData = readPairVectorFromFile(DQTpath);
    list<std::pair<double, double>> betasAndDiffs;
    for (int i = 0; i < EDData.size(); i++) {
        double absDiff = EDData[i].second - DQTData[i].second;
        double relDiff;
        relDiff = absDiff / EDData[i].second;
        /*if (EDData[i].second > 1e-10 && absDiff > 1e-10) {
            relDiff = absDiff / EDData[i].second;
        } else {
            relDiff = 0.0;
        }*/
        betasAndDiffs.emplace_back(std::pair<double, double>( EDData[i].first, abs( relDiff)));
    }
    savePairsToFile(betasAndDiffs, outPath);
}

// Reads double pairs from a file.
vector<std::pair<double, double>> readPairVectorFromFile(const std::string & path) {
    vector<std::pair<double, double>> vals;
    std::ifstream file;
    file.open(path);
    std::string line;
    while (getline(file, line)) {
        std::vector<double> tempV;
        std::istringstream lineStream(line);
        string temp;
        while (std::getline(lineStream, temp, ' ')) {
            tempV.emplace_back( stod(temp) );
        }
        vals.emplace_back( std::pair<double, double>(tempV[0], tempV[1]) );
    }

    return vals;
}

// Normalises a list of complex vectors to be interpreted as one large vector.
void normaliseListOfVectors(vector<VectorXcd> & vec) {
    double norm2 = 0;
    for (VectorXcd & v : vec) {
        norm2 += pow(v.norm(),2);
    }
    double norm = sqrt(norm2);
    for (VectorXcd & v : vec) {
        v /= norm;
    }
}

// Generates a random complex vector using a true random seed for a pseudo random number generator.
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

// Generates a random complex vector using a true random seed for a pseudo random number generator.
VectorXcd randomComplexVector(int vecSize, double stdDev) {
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

    return psi;
}
