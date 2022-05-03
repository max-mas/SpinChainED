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
    double j_ratio = 2;


    MatrixXd H = naiveHamiltonian(j_ratio, N);
    //printMatrix(H);
    Eigen::VectorXd erg = H.eigenvalues().real();
    std::sort(erg.begin(), erg.end());
    printEnergies(erg);
    //printMatrix(H);

    std::cout << std::endl;

    list<MatrixXd> H2 = magnetizationHamiltonian(j_ratio, N);
    vector<double> erg2 = getEnergiesFromBlocks(H2, N);
    printEnergies(erg2);

    std::cout << std::endl;

    list<list<MatrixXcd>> H3 = momentumHamiltonian(j_ratio, N);
    vector<double> erg3 = getEnergiesFromBlocks(H3, N);
    printEnergies(erg3);
    //list<MatrixXcd> H3_1 = blkdiag(H3);
    //MatrixXcd H3_2 = blkdiag(H3_1, pow(2, N));
    //printMatrix(H3_2);
    //MatrixXcd H3_2_tr = H3_2.transpose();
    //printMatrix(H3_2_tr);

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
        for (int i = 0; i < N; i++) {
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
            for (int i = 0; i < N; i++) {
                int j = (i+2) % N;
                if (getBit(a, i) == getBit(a, j)) {
                    H_subspace_list.back()(k, k) += J_ratio* 0.25;
                } else {
                    H_subspace_list.back()(k, k) += J_ratio * -0.25;

                    int b = a;
                    flipBit(b, i);
                    flipBit(b, j);
                    int l = findState(s_vector_m, b);
                    H_subspace_list.back()(k, l) = J_ratio * 0.5;
                }
            }
        }
    }
    return H_subspace_list;
}

list<list<MatrixXcd>> momentumHamiltonian(double J_ratio, int N) {
    list<list<MatrixXcd>> H_subspace_list;

    //int H_index = 0;
    for (int m_setter = 0; m_setter <= N; m_setter++) {
        double m = -N/2.0 + m_setter;
        int n_up = round( m + N/2.0);

        vector<int> s_vector_m;
        for (int s = 0; s < pow(2, N); s++) {
            if (bitSum(s) == n_up) { s_vector_m.push_back(s); }
        }
        int M = s_vector_m.size();

        list<MatrixXcd> H_subSubspace_list(N);

        for (int k = -trunc((N+2)/4) + 1 ; k <= trunc(N/4); k++ ) {
            vector<int> s_vector_k, R_vector;
            for (int i = 0; i < M; i++) {
                int s = s_vector_m[i];
                int R = checkState(s, k, N);
                if (R >= 0) {
                    s_vector_k.push_back(s);
                    R_vector.push_back(R);
                }
            }
            int K = s_vector_k.size();

            H_subSubspace_list.emplace_back(MatrixXcd::Zero(K, K));

            for (int l = 0; l < K; l++) {
                int a = s_vector_k[l];
                for (int i = 0; i < N; i++) {
                    int j = (i+1) % N;
                    if (getBit(a, i) == getBit(a, j)) {
                        H_subSubspace_list.back()(l, l) += 0.25;
                    } else {
                        H_subSubspace_list.back()(l, l) += -0.25;

                        int b = a;
                        flipBit(b, i);
                        flipBit(b, j);
                        vector<int> r_L = representative(b, N);
                        int f = findState(s_vector_k, r_L[0]);

                        if (f >= 0) {
                            complex<double> offDiagEl = 0.5 * std::sqrt((double)R_vector[l] / (double) R_vector[f])
                                    * std::exp(complex<double>(0,4.0 * M_PI * (double) k * (double) r_L[1] / (double) N));
                            H_subSubspace_list.back()(l, f) += offDiagEl;
                        }
                    }
                    int n = (i+2) % N;
                    if (getBit(a, i) == getBit(a, n) && i < N) {
                        H_subSubspace_list.back()(l, l) += J_ratio * 0.25;
                    } else if (i < N) {
                        H_subSubspace_list.back()(l, l) += J_ratio * -0.25;

                        int b = a;
                        flipBit(b, i);
                        flipBit(b, n);
                        vector<int> r_L = representative(b, N);
                        int f = findState(s_vector_k, r_L[0]);

                        if (f >= 0) {
                            complex<double> offDiagEl = 0.5 * std::sqrt((double)R_vector[l] / (double) R_vector[f])
                                    * std::exp( complex<double>(0,4.0 * M_PI * (double) k * (double) r_L[1] / (double) N));
                            H_subSubspace_list.back()(l, f) += J_ratio * offDiagEl;
                        }
                    }
                }
            }
        }
        H_subspace_list.push_back(H_subSubspace_list);
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


MatrixXcd blkdiag(const list<MatrixXcd> & matrix_list, int totalSize) {
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

list<MatrixXcd> blkdiag(const list<list<MatrixXcd>> & matrix_doubleList) {
    list<MatrixXcd> bdmlist;
    for (list<MatrixXcd> l : matrix_doubleList) {
        int size = 0;
        for (MatrixXcd m : l) {
            size += m.rows();
        }
        bdmlist.emplace_back( blkdiag(l, size) );
    }
    return bdmlist;
}