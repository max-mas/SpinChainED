#include "momentumHamiltonian_sparse.h"

using std::vector;
using std::complex;

using Eigen::SparseMatrix;

SparseMatrix<complex<double>> momentumHamiltonian_sparse(double J_ratio, int N) {

    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    typedef Eigen::Triplet<complex<double>> matFil;

    // init list of elements to be filled
    vector<matFil> elementList;

    // "global" index corresponding to the (1,1)-Element of the current block
    int currentMatRootIndex = 0;

    // loop over all magnetizations m
    for (int m_setter = 0; m_setter <= N; m_setter++) {
        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, m_setter);
        int M = s_vector_m.size();

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

            // fill elements of (m,k)-block

            for (int l = 0; l < K; l++) {
                int a = s_vector_k[l];
                //std::cout << a << " " << k << " " << R_vector[l] << std::endl;
                for (int i = 0; i < N; i++) {
                    int j = (i+1) % N;
                    if (getBit(a, i) == getBit(a, j)) {
                        elementList.emplace_back(matFil(l + currentMatRootIndex, l + currentMatRootIndex, 0.25));
                    } else {
                        elementList.emplace_back(matFil(l + currentMatRootIndex, l + currentMatRootIndex, -0.25));

                        int b = a;
                        flipBit(b, i);
                        flipBit(b, j);
                        vector<int> r_L = representative(b, N);
                        int f = findState(s_vector_k, r_L[0]);

                        if (f >= 0) {
                            complex<double> offDiagEl = 0.5 * sqrt((double)R_vector[l] / (double) R_vector[f])
                                                        * std::exp(complex<double>(0,4.0 * M_PI * (double) k * (double) r_L[1] / (double) N));
                            elementList.emplace_back(matFil(l + currentMatRootIndex, f + currentMatRootIndex, offDiagEl));
                        }
                    }
                    int n = (i+2) % N;
                    if (getBit(a, i) == getBit(a, n) && i < N) {
                        elementList.emplace_back(matFil(l + currentMatRootIndex, l + currentMatRootIndex, J_ratio * 0.25));
                    } else if (i < N) {
                        elementList.emplace_back(matFil(l + currentMatRootIndex, l + currentMatRootIndex, J_ratio * -0.25));

                        int b = a;
                        flipBit(b, i);
                        flipBit(b, n);
                        vector<int> r_L = representative(b, N);
                        int f = findState(s_vector_k, r_L[0]);

                        if (f >= 0) {
                            complex<double> offDiagEl = 0.5 * sqrt((double)R_vector[l] / (double) R_vector[f])
                                                        * std::exp( complex<double>(0,4.0 * M_PI * (double) k * (double) r_L[1] / (double) N));
                            elementList.emplace_back(matFil(l + currentMatRootIndex, f + currentMatRootIndex, J_ratio * offDiagEl));
                        }
                    }
                }
            }

            currentMatRootIndex += K;
        }

    }

    SparseMatrix<complex<double>> H(pow(2, N), pow(2, N));
    H.setFromTriplets(elementList.begin(), elementList.end());
    H.makeCompressed();

    return H;
}
