//
// Created by mmaschke on 03/05/22.
//

#include "hamiltonianBuilders.h"

using std::vector;
using std::complex;
using std::list;
using std::pow;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;

MatrixXd naiveHamiltonian(double J_ratio, int N) {
    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    // init Matrix
    int matrixSize = pow(2, N);
    MatrixXd H = MatrixXd::Zero(matrixSize, matrixSize);

    // set Matrix Elements for all states
    for (int a = 0; a < matrixSize; a++) {
        setHElement_naive(J_ratio, N, H, a);
    }
    return H;
}

void setHElement_naive(double J_ratio, int N, MatrixXd & H, int a) {
    // J2 coupling
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
    // J1 coupling
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

list<MatrixXd> magnetizationHamiltonian(double J_ratio, int N) {
    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    // init empty block list
    list<MatrixXd> H_subspace_list;

    // loop over all magnetizations
    for (int m_setter = 0; m_setter <= N; m_setter++) {
        // Calculate magnetization and number of "up"-states for given magnetization.
        double m = -N/2.0 + m_setter;
        int n_up = round( m + N/2.0);

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, n_up);
        int M = s_vector_m.size();

        // Generate Block with correct size and fill elements.
        H_subspace_list.emplace_back(MatrixXd::Zero(M,M));
        for (int k = 0; k < M; k++) {
            setHElement_magnetization(J_ratio, N, H_subspace_list, s_vector_m, k);
        }
    }
    return H_subspace_list;
}

vector<int> getStates_m(int N, int n_up) {
    vector<int> s_vector_m;
    for (int s = 0; s < pow(2, N); s++) {
        if (bitSum(s) == n_up) {s_vector_m.push_back(s);}
    }
    return s_vector_m;
}

void setHElement_magnetization(double J_ratio, int N, list<MatrixXd> &H_subspace_list, const vector<int> &s_vector_m,
                               int k) {
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

list<list<MatrixXcd>> momentumHamiltonian(double J_ratio, int N) {
    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    // init list of block lists
    list<list<MatrixXcd>> H_subspace_list;

    // loop over all magnetizations m
    for (int m_setter = 0; m_setter <= N; m_setter++) {
        // Calculate magnetization and number of "up"-states for given magnetization.
        double m = -N/2.0 + m_setter;
        int n_up = round( m + N/2.0);

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, n_up);
        int M = s_vector_m.size();

        // init list of blocks
        list<MatrixXcd> H_subSubspace_list(N);

        // loop over all possible momenta k
        for (int k = -trunc((N+2)/4) + 1 ; k <= trunc(N/4); k++ ) {

            // find k-compatible states and their periodicities
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

            // init (m,k)-block
            H_subSubspace_list.emplace_back(MatrixXcd::Zero(K, K));

            // fill elements of (m,k)-block
            setHElement_momentum(J_ratio, N, H_subSubspace_list, k, s_vector_k, R_vector, K);
        }
        H_subspace_list.push_back(H_subSubspace_list);
    }
    return H_subspace_list;
}

void setHElement_momentum(double J_ratio, int N, list<MatrixXcd> &H_subSubspace_list, int k,
                          const vector<int> &s_vector_k, const vector<int> &R_vector, int K) {
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
                    complex<double> offDiagEl = 0.5 * sqrt((double)R_vector[l] / (double) R_vector[f])
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
                    complex<double> offDiagEl = 0.5 * sqrt((double)R_vector[l] / (double) R_vector[f])
                                                * std::exp( complex<double>(0,4.0 * M_PI * (double) k * (double) r_L[1] / (double) N));
                    H_subSubspace_list.back()(l, f) += J_ratio * offDiagEl;
                }
            }
        }
    }
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
