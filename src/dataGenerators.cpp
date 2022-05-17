/**
 * This cpp file contains functions used for generating data of interest (excitation energies, specific heat and
 * susceptibilities) for later plotting in Python.
*/

#include "dataGenerators.h"

#include <utility>

using std::complex;
using std::list;
using std::pow;
using std::vector;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;

// Saves excitation energy (i.e. energy difference between ground and 1st excited state) vor varying values of J1/J2.
// Note: If ground state is degenerate, zero is returned.
void saveExcitationErgsForVaryingJ(int N, int dataPointNum, double start, double end, const std::string & path) {
    Eigen::VectorXd J_ratios = Eigen::VectorXd::LinSpaced(dataPointNum, start, end);
    vector<double> diffs;
#pragma omp parallel for default(none) shared(diffs, J_ratios, N)
    for (int i = 0; i < J_ratios.size(); i++) {
        if (N % 4 == 0 && N >= 8) {
            list<list<list<MatrixXd>>> H = spinInversionHamiltonian(J_ratios[i], N, 0, N);
            vector<double> erg = getEnergiesFromBlocks(H, true);
            writeThreadSafe(diffs, {abs(erg[0]-erg[1])});
        } else if (N % 2 == 0 && N >= 6) {
            list<list<MatrixXcd>> H = momentumHamiltonian(J_ratios[i], N);
            vector<double> erg = getEnergiesFromBlocks(H, true);
            writeThreadSafe(diffs, {abs(erg[0]-erg[1])});
        }
    }
    list<std::pair<double, double>> out;
    for (int i = 0; i < diffs.size(); i++) {
        out.emplace_back( std::pair<double, double>(J_ratios[i], diffs[i]) );
    }
    savePairsToFile(out, path);
}

// Saves specific heat at given temperature/beta for varying values of J1/J2.
void saveSpecificHeatForVaryingJ(int N, int dataPointNum, double betaOrT, double start, double end, bool isBeta, std::string path) {
    Eigen::VectorXd J_ratios = Eigen::VectorXd::LinSpaced(dataPointNum, start, end);
    vector<double> C;
    for (double J_ratio : J_ratios) {
        if (N % 4 == 0 && N >= 8) {
            list<list<list<MatrixXd>>> H = spinInversionHamiltonian(J_ratio, N, 0, N);
            vector<double> erg = getParityErgsThreaded(H, N, true);
            C.emplace_back(specificHeat(erg, betaOrT, isBeta) / N);
        } else if (N % 2 == 0 && N >= 6) {
            list<list<MatrixXcd>> H = momentumHamiltonian(J_ratio, N);
            vector<double> erg = getMomentumErgsThreaded(H, N, true);
            C.emplace_back(specificHeat(erg, betaOrT, isBeta) / N);
        }
    }
    list<std::pair<double, double>> out;
    for (int i = 0; i < C.size(); i++) {
        out.emplace_back( std::pair<double, double>(J_ratios[i], C[i]) );
    }
    savePairsToFile(out, path);
}

// Saves specific heat for a given value of J1/J2 for varying temperature/beta.
void saveSpecificHeatsForVaryingTemp(int N, int dataPointNum, double J_ratio, double start, double end,
                                     bool isBeta, std::string path) {
    Eigen::VectorXd Ts = Eigen::VectorXd::LinSpaced(dataPointNum, start, end);
    vector<double> C;

    if (N % 4 == 0 && N >= 8) {
        list<list<list<MatrixXd>>> H = spinInversionHamiltonian(J_ratio, N, 0, N);
        vector<double> erg = getParityErgsThreaded(H, N, true);
        for (int i = 0; i < Ts.size(); i++) {
            C.emplace_back(specificHeat(erg, Ts[i], isBeta) / N);
        }
    } else if (N % 2 == 0 && N >= 6) {
        list<list<MatrixXcd>> H = momentumHamiltonian(J_ratio, N);
        vector<double> erg = getMomentumErgsThreaded(H, N, true);
        for (int i = 0; i < Ts.size(); i++) {
            C.emplace_back(specificHeat(erg, Ts[i], isBeta) / N);
        }
    }

    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(Ts[j], C[j]) );
    }
    savePairsToFile(out, path);
}

// Saves susceptibilities at a given temperature/beta for varying values of J1/J2.
void saveSusceptibilitiesForVaryingJ(int N, int dataPointNum, double betaOrT, double start, double end,
                                     bool isBeta, std::string path) {
    Eigen::VectorXd J_ratios = Eigen::VectorXd::LinSpaced(dataPointNum, start, end);
    vector<int> states = getStates_m(N, N/2);
    MatrixXd S_2 = spinOperator_sq(states, N);

    vector<double> susceptibilities;
    for (double J_ratio : J_ratios) {
        if (N % 4 == 0 && N >= 8) {
            list<list<list<MatrixXd>>> H_m0 = spinInversionHamiltonian(J_ratio,  N, N/2, N/2);
            vector<double> ergs;
            MatrixXd U = buildTransformMatrix_parity(H_m0, ergs);

            MatrixXd T = U.transpose() * S_2 * U;

            susceptibilities.emplace_back( susceptibility(ergs, betaOrT, isBeta, T) / (double) N );
        } else if (N % 2 == 0 && N >= 6) {
            MatrixXd H_m0 = getMagnetizationBlock(J_ratio, 0, N);
            Eigen::SelfAdjointEigenSolver<MatrixXd> sol(H_m0);
            Eigen::VectorXd erg = sol.eigenvalues().real();
            vector<double> ergs_stl(erg.data(), erg.data() + erg.size());
            const Eigen::MatrixXcd & U = sol.eigenvectors();

            MatrixXcd T = U.adjoint() * S_2 * U;

            susceptibilities.emplace_back( susceptibility(ergs_stl, betaOrT, isBeta, T) / (double) N );
        }
    }
    list<std::pair<double, double>> out;
    for (int i = 0; i < susceptibilities.size(); i++) {
        out.emplace_back( std::pair<double, double>(J_ratios[i], susceptibilities[i]) );
    }
    savePairsToFile(out, path);
}

// Saves susceptibilities for a given value of J1/J2 for varying temperature/beta.
void saveSusceptibilitesForVaryingTemp(int N, int dataPointNum, double J_ratio, double start, double end,
                                       bool isBeta, std::string path) {
    Eigen::VectorXd Ts = Eigen::VectorXd::LinSpaced(dataPointNum, start, end);
    vector<int> states = getStates_m(N, N/2);
    MatrixXd S_2 = spinOperator_sq(states, N);

    vector<double> susceptibilities(dataPointNum);
    if (N % 4 == 0 && N >= 8) {
        list<list<list<MatrixXd>>> H_m0 = spinInversionHamiltonian(J_ratio,  N, N/2, N/2);
        vector<double> ergs;
        MatrixXd U = buildTransformMatrix_parity(H_m0, ergs);

        MatrixXd T = U.transpose() * S_2 * U;

        for (int i = 0; i < dataPointNum; i++) {
            susceptibilities[i] =  susceptibility(ergs, Ts[i], isBeta, T) / N ;
        }
    } else if (N % 2 == 0 && N >= 6) {
        MatrixXd H_m0 = getMagnetizationBlock(J_ratio, 0, N);
        Eigen::SelfAdjointEigenSolver<MatrixXd> sol(H_m0);
        Eigen::VectorXd erg = sol.eigenvalues().real();
        vector<double> ergs_stl(erg.data(), erg.data() + erg.size());
        const Eigen::MatrixXcd & U = sol.eigenvectors();

        MatrixXcd T = U.adjoint() * S_2 * U;

        for (int i = 0; i < dataPointNum; i++) {
            susceptibilities[i] =  susceptibility(ergs_stl, Ts[i], isBeta, T) / N ;
        }
    }

    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(Ts[j], susceptibilities[j]) );
    }
    savePairsToFile(out, path);
}

// Saves tuples of (E, k) to file.
void saveEnergyDispersion(int N, double J_ratio, std::string path) {
    list<list<MatrixXcd>> H_list = momentumHamiltonian(J_ratio, N);
    vector<vector<vector<double>>> ergs = getEnergiesFromBlocksByK(H_list);
    list<std::pair<int, double>> out;
    for (vector<vector<double>> m_block : ergs) {
        int k = -trunc((N+2)/4) + 1;
        for (vector<double> k_block : m_block) {
            for (double erg : k_block) {
                out.emplace_back( std::pair<int, double>(k, erg) );
            }
            k++;
        }
    }
    savePairsToFile(out, path);
}

// Saves tuples of (E, m, k) to file.
void saveEnergyDispersionWithMag(int N, double J_ratio, std::string path) {
    list<list<MatrixXcd>> H_list = momentumHamiltonian(J_ratio, N);
    vector<vector<vector<double>>> ergs = getEnergiesFromBlocksByK(H_list);
    list<std::tuple<int, int, double>> out;
    int m = -N/2.0;
    for (vector<vector<double>> m_block : ergs) {
        int k = -trunc((N+2)/4) + 1;
        for (vector<double> k_block : m_block) {
            for (double erg : k_block) {
                out.emplace_back( std::tuple<int, int, double>(m, k, erg) );
            }
            k++;
        }
        m++;
    }
    saveTripleToFile(out, path);
}

// Used to write data tuples to a file at path.
template <typename T, typename U>
void savePairsToFile(const list<std::pair<T, U>> & pairList, std::string path) {
    std::ofstream File;
    File.open(path);
    for (std::pair<T, U> p : pairList) {
        File << p.first << " " << p.second << "\n";
    }
    File.close();
}

template <typename T, typename U, typename V>
void saveTripleToFile(const list<std::tuple<T, U, V>> & pairList, std::string path) {
    std::ofstream File;
    File.open(path);
    for (std::tuple<T, U, V> p : pairList) {
        File << std::get<0>(p) << " " << std::get<1>(p) << " " << std::get<2>(p) << "\n";
    }
    File.close();
}

// Enables threaded writing to a vector.
void writeThreadSafe (vector<vector<double>> & writeTo, const vector<double> & writeFrom) {
    static std::mutex mu;
    std::lock_guard<std::mutex> lock(mu);
    writeTo.emplace_back(writeFrom);
}

// Enables threaded writing to a vector.
void writeThreadSafe (vector<double> & writeTo, const vector<double> & writeFrom) {
    static std::mutex mu;
    std::lock_guard<std::mutex> lock(mu);
    for (double d : writeFrom) {
        writeTo.emplace_back(d);
    }
}

MatrixXd buildTransformMatrix_parity(const list<list<list<MatrixXd>>> & H_m0, vector<double> & ergs) {
    int totalSize = 0;
    list<list<MatrixXd>> EV_blocks;
    for (const list<list<MatrixXd>> & l : H_m0) { // if this ever has more than one element, kaboom
        for (const list<MatrixXd> & m : l) {
            list<MatrixXd> EV_subBlocks;
            for (const MatrixXd & mat : m) {
                totalSize += mat.cols();

                Eigen::SelfAdjointEigenSolver<MatrixXd> sol(mat);
                Eigen::VectorXd blockErgs = sol.eigenvalues().real();
                vector<double> ergs_stl(blockErgs.begin(), blockErgs.end());
                writeThreadSafe(ergs, ergs_stl);

                MatrixXd U_block = sol.eigenvectors();
                EV_subBlocks.emplace_back(U_block);
            }
            EV_blocks.emplace_back(EV_subBlocks);
        }
    }

    list<MatrixXd> temp = blkdiag(EV_blocks);
    MatrixXd ret = blkdiag(temp, totalSize);

    return ret;
}