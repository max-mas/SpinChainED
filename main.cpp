#include "main.h"
#include <chrono>

using std::vector;
using std::complex;
using std::list;
using std::pow;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;

int main(int argc, char* argv[]) {

    int N = 6;

    MatrixXd H = naiveHamiltonian(1, N);
    //printMatrix(H);
    Eigen::VectorXd erg = H.eigenvalues().real();
    std::sort(erg.begin(), erg.end());
    printEnergies(erg);

    std::cout << std::endl;

    list<MatrixXd> H2 = magnetizationHamiltonian(1, N);
    vector<double> erg2 = getEnergiesFromBlocks(H2, N);
    printEnergies(erg2);

    return 0;
}

MatrixXd naiveHamiltonian(double J_ratio, int N) {
    int matrixSize = pow(2, N); // N must be even or this no longer describes the correct system.
    MatrixXd H = MatrixXd::Zero(matrixSize, matrixSize);

    for (int a = 0; a < matrixSize; a++) {
        for (int i = 0; i < N; i++) {
            int j = (i+1) % N;
            int k = (i+2) % N;
            if (getBit(a, i) == getBit(a, j)) {
                H(a, a) += 0.25;
            } else {
                H(a, a) += -0.25;
                int b = a;
                flipBit(b, i);
                flipBit(b, j);
                H(a, b) = 0.5;
            }
        }
        for (int i = 0; i < N-1; i++) {
            int j = (i+2) % N;
            if (getBit(a, i) == getBit(a, j)) {
                H(a, a) += J_ratio *  0.25;
            } else {
                H(a, a) += J_ratio * -0.25;
                int b = a;
                flipBit(b, i);
                flipBit(b, j);
                H(a, b) = J_ratio * 0.5;
            }
        }
    }
    return H;
}

list<MatrixXd> magnetizationHamiltonian(double J_ratio, int N) {
    list<MatrixXd> H_subspace_list;

    for (int m_setter = 0; m_setter <= N; m_setter++) {
        // Calculate magnetization and number of "up"-states for given magnetization.
        double m = -N/2.0 + m_setter;
        int n_up = round( m + N/2.0);

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m;
        for (int s = 0; s < pow(2, N); s++) {
            if (bitSum(s) == n_up) {s_vector_m.push_back(s);}
        }
        int M = s_vector_m.size();

        // Generate Block with correct size and fill elements.
        H_subspace_list.emplace_back(MatrixXd::Zero(M,M));
        for (int k = 0; k < M; k++) {
            int a = s_vector_m[k];
            for (int i = 0; i <N; i++) {
                int j = (i+1) % N;
                if (getBit(a, i) == getBit(a, j)) {
                    H_subspace_list.back()(k, k) += 0.25;
                } else {
                    H_subspace_list.back()(k, k) += -0.25;

                    int b = a;
                    flipBit(b, i);
                    flipBit(b, j);
                    int l = findState(s_vector_m, b);
                    H_subspace_list.back()(k, l) = 0.5;
                }
            }
            for (int i = 0; i <N-1; i++) {
                int j = (i+2) % N;
                if (getBit(a, i) == getBit(a, j)) {
                    H_subspace_list.back()(k, k) += 0.25;
                } else {
                    H_subspace_list.back()(k, k) += -0.25;

                    int b = a;
                    flipBit(b, i);
                    flipBit(b, j);
                    int l = findState(s_vector_m, b);
                    H_subspace_list.back()(k, l) = 0.5;
                }
            }
        }
    }
    return H_subspace_list;
}

vector<double> getEnergiesFromBlocks(const list<MatrixXcd> & H_list, int N) {
    vector<double> energies(pow(2, N), 0);
    int j = 0;
    for (const MatrixXcd & mat : H_list) {
        Eigen::VectorXcd blockEnergies = mat.eigenvalues();
        std::for_each(blockEnergies.begin(), blockEnergies.end(),
                      [&](complex<double> & d){energies[j] = d.real(); j++;});
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<MatrixXd> & H_list, int N) {
    vector<double> energies(pow(2, N), 0);
    int j = 0;
    for (const MatrixXd & mat : H_list) {
        Eigen::VectorXcd blockEnergies = mat.eigenvalues();
        std::for_each(blockEnergies.begin(), blockEnergies.end(),
                      [&](complex<double> & d){energies[j] = d.real(); j++;});
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<list<MatrixXcd>> & H_list, int N) {
    vector<double> energies(pow(2, N), 0);
    int j = 0;
    for (const list<MatrixXcd> & subList : H_list) {
        for (const MatrixXcd & mat : subList) {
            Eigen::VectorXcd blockEnergies = mat.eigenvalues();
            std::for_each(blockEnergies.begin(), blockEnergies.end(),
                          [&](complex<double> &d) {
                              energies[j] = d.real();
                              j++; });
        }
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}


MatrixXcd blkdiag(const vector<MatrixXcd> & matrix_list, int totalSize) {
    MatrixXcd bdm = MatrixXcd::Zero(totalSize, totalSize);
    int curr_index = 0;
    for (const MatrixXcd & mat : matrix_list)
    {
        if (mat.cols() == 0) {
            continue;
        }
        bdm.block(curr_index, curr_index, mat.rows(), mat.cols()) = mat;
        curr_index += (int) mat.rows();
    }
    return bdm;
}

MatrixXd blkdiag(const vector<MatrixXd> & matrix_list, int totalSize) {
    MatrixXd bdm = MatrixXd::Zero(totalSize, totalSize);
    int curr_index = 0;
    for (const MatrixXd & mat : matrix_list)
    {
        if (mat.cols() == 0) {
            continue;
        }
        bdm.block(curr_index, curr_index, mat.rows(), mat.cols()) = mat;
        curr_index += (int) mat.rows();
    }
    return bdm;
}




