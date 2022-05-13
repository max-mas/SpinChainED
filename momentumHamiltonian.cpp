//
// Created by MaxM on 13/05/2022.
//

#include "momentumHamiltonian.h"

using std::vector;
using std::list;
using std::complex;

using Eigen::MatrixXcd;

// Generate list of list of (m, k)-blocks of the Hamiltonian for a given N.
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

// Set matrix element for a given state (momentum-approach).
void setHElement_momentum(double J_ratio, int N, list<MatrixXcd> &H_subSubspace_list, int k,
                          const vector<int> &s_vector_k, const vector<int> &R_vector, int K) {
    for (int l = 0; l < K; l++) {
        int a = s_vector_k[l];
        //std::cout << a << " " << k << " " << R_vector[l] << std::endl;
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

// Checks if state int s is smallest among its translations and compatible with momentum k.
int checkState(const int s, const int k, const int N) {
    int t = s;
    for (int i = 1; i <= N/2; i++) {
        cycleBits2(t, N); //translate state
        if (t < s) {return -1;}
        else if (t == s) {
            if ( k % (int) std::trunc((double )N/(2*i)) ) {return -1;} // check compatibility with k
            return i;
        }
    }
    return -1;
}


// Generates representative for given state int s. 1st return: representative. 2nd return:
// Needed number of translations.
std::vector<int> representative(const int s, const int N) {
    int r = s;
    int t = s;
    int l = 0;

    for (int i = 1; i < N/2; i++) {
        cycleBits2(t, N);
        if (t < r) {r = t; l = i;}
    }
    return {r, l};
}
