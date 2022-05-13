//
// Created by MaxM on 13/05/2022.
//

#include "magnetizationHamiltonian.h"

using Eigen::MatrixXd;
using std::vector;
using std::list;

// For a given m, get the corresponding block of the Hamiltonian.
MatrixXd getMagnetizationBlock(double J_ratio, double m, int N) {
    int n_up = round( m + N/2.0);
    vector<int> s_vector_m = getStates_m(N, n_up);
    int M = s_vector_m.size();

    MatrixXd m_block =MatrixXd::Zero(M,M);
    for (int k = 0; k < M; k++) {
        setHElement_magnetization(J_ratio, N, m_block, s_vector_m, k);
    }
    return m_block;
}

// Get a list of all magnetization-blocks for a given N.
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
            setHElement_magnetization(J_ratio, N, H_subspace_list.back(), s_vector_m, k);
        }
    }
    return H_subspace_list;
}

// Set matrix element for a given state (magnetization approach).
void setHElement_magnetization(double J_ratio, int N, MatrixXd & H_subspace_block, const vector<int> &s_vector_m,
                               int k) {
    int a = s_vector_m[k];
    for (int i = 0; i <N; i++) {
        int j = (i+1) % N;
        if (getBit(a, i) == getBit(a, j)) {
            H_subspace_block(k, k) += 0.25;
        } else {
            H_subspace_block(k, k) += -0.25;

            int b = a;
            flipBit(b, i);
            flipBit(b, j);
            int l = findState(s_vector_m, b);
            H_subspace_block(k, l) = 0.5;
        }
    }
    for (int i = 0; i < N; i++) {
        int j = (i+2) % N;
        if (getBit(a, i) == getBit(a, j)) {
            H_subspace_block(k, k) += J_ratio* 0.25;
        } else {
            H_subspace_block(k, k) += J_ratio * -0.25;

            int b = a;
            flipBit(b, i);
            flipBit(b, j);
            int l = findState(s_vector_m, b);
            H_subspace_block(k, l) = J_ratio * 0.5;
        }
    }
}
