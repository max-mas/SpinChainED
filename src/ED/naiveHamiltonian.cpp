/**
* This file contains methods used to construct the full-sized Hamiltonian as well as the total-spin operator S².
*/

#include "naiveHamiltonian.h"

using std::vector;
using std::complex;
using std::list;
using std::pow;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;

// Generate full-size S² for given N.
MatrixXd spinOperator_sq(int N) {
    int matSize = pow(2, N);
    MatrixXd S_2 = N*0.75*MatrixXd::Identity(matSize, matSize);
    for (int a = 0; a < matSize; a++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                if (getBit(a, i) == getBit(a, j)) {
                    S_2(a, a) += 0.5;
                } else {
                    S_2(a, a) += -0.5;
                    int b = a;
                    flipBit(b, i);
                    flipBit(b, j);
                    S_2(a, b) += 1;
                }
            }
        }
    }
    return S_2;
}

// Generate S² for a given subset of states (in practice: m = 0 states).
MatrixXd spinOperator_sq(vector<int> states, int N) {
    int matSize = states.size();
    MatrixXd S_2 = N*0.75*MatrixXd::Identity(matSize, matSize);
    for (int k = 0; k < matSize; k++) {
        int a = states[k];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                if (getBit(a, i) == getBit(a, j)) {
                    S_2(k, k) += 0.5;
                } else {
                    S_2(k, k) += -0.5;

                    int b = a;
                    flipBit(b, i);
                    flipBit(b, j);
                    int l = findState(states, b);
                    S_2(k, l) += 1;
                }
            }
        }
    }
    return S_2;
}

// Generate full-size hamiltonian using the naive approach.
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

// Set matrix element for a given state (naive approach).
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