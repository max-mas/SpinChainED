/**
 * This cpp file contains functions used for generating data of interest (excitation energies, specific heat and
 * susceptibilities) for later plotting in Python.
*/

#include "dataGenerators.h"

using std::complex;
using std::list;
using std::pow;
using std::vector;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;

// Saves excitation energy (i.e. energy difference between ground and 1st excited state) vor varying values of J1/J2.
// Note: If ground state is degenerate, zero is returned.
void saveExcitationErgsForVaryingJ(int N, int dataPointNum, double start, double end, std::string path) {
    Eigen::VectorXd J_ratios = Eigen::VectorXd::LinSpaced(dataPointNum, start, end);
    vector<double> diffs;
    for (double J_ratio : J_ratios) {
        list<list<MatrixXcd>> H = momentumHamiltonian(J_ratio, N);
        vector<double> erg = getEnergiesFromBlocks(H, N);
        diffs.emplace_back( abs(erg[0]-erg[1]) );
    }
    list<std::pair<double, double>> out;
    for (int i = 0; i < diffs.size(); i++) {
        out.emplace_back( std::pair<double, double>(J_ratios[i], diffs[i]) );
    }
    savePairsToFile(out, path);
}

void saveExcitationEnergiesByK() {

}

// Saves specific heat at given temperature/beta for varying values of J1/J2.
void saveSpecificHeatForVaryingJ(int N, int dataPointNum, double betaOrT, double start, double end, bool isBeta, std::string path) {
    Eigen::VectorXd J_ratios = Eigen::VectorXd::LinSpaced(dataPointNum, start, end);
    vector<double> C;
    for (double J_ratio : J_ratios) {
        list<list<MatrixXcd>> H = momentumHamiltonian(J_ratio, N);
        vector<double> erg = getMomentumErgsThreaded(H, N);
        C.emplace_back(specificHeat(erg, betaOrT, isBeta) / N);
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
    list<list<MatrixXcd>> H = momentumHamiltonian(J_ratio, N);
    vector<double> erg = getMomentumErgsThreaded(H, N);
    for (int i = 0; i < Ts.size(); i++) {
        C.emplace_back( specificHeat(erg, Ts[i], isBeta) /N );
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

        MatrixXd H_m0 = getMagnetizationBlock(J_ratio, 0, N);
        Eigen::SelfAdjointEigenSolver<MatrixXd> sol(H_m0);
        Eigen::VectorXd erg = sol.eigenvalues().real();
        vector<double> ergs_stl(erg.data(), erg.data() + erg.size());
        const Eigen::MatrixXcd & U = sol.eigenvectors();

        MatrixXcd T = U.adjoint() * H_m0 * U;

        susceptibilities.emplace_back( susceptibility(ergs_stl, betaOrT, isBeta, T) / (double) N );
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

    MatrixXd H_m0 = getMagnetizationBlock(J_ratio, 0, N);
    Eigen::SelfAdjointEigenSolver<MatrixXd> sol(H_m0);
    Eigen::VectorXd erg = sol.eigenvalues().real();
    vector<double> ergs_stl(erg.data(), erg.data() + erg.size());
    const Eigen::MatrixXcd & U = sol.eigenvectors();

    MatrixXcd T = U.adjoint() * H_m0 * U;

    vector<double> susceptibilities(dataPointNum);
    for (int i = 0; i < dataPointNum; i++) {
        susceptibilities[i] =  susceptibility(ergs_stl, Ts[i], isBeta, T) / N ;
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

