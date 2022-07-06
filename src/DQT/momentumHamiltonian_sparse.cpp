#include "momentumHamiltonian_sparse.h"

using std::vector;
using std::list;
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

SparseMatrix<complex<double>> spinOp2_momentum_sparse(int N) {

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
                    for (int j = 0; j < i; j++) {
                        if (getBit(a, i) == getBit(a, j)) {
                            elementList.emplace_back(matFil(l + currentMatRootIndex, l + currentMatRootIndex, 0.5));
                        } else {
                            elementList.emplace_back(matFil(l + currentMatRootIndex, l + currentMatRootIndex, -0.5));

                            int b = a;
                            flipBit(b, i);
                            flipBit(b, j);
                            vector<int> r_L = representative(b, N);
                            int f = findState(s_vector_k, r_L[0]);

                            if (f >= 0) {
                                complex<double> offDiagEl = sqrt((double)R_vector[l] / (double) R_vector[f])
                                                            * std::exp(complex<double>(0,4.0 * M_PI * (double) k * (double) r_L[1] / (double) N));
                                elementList.emplace_back(matFil(l + currentMatRootIndex, f + currentMatRootIndex, offDiagEl));
                            }
                        }
                    }
                }
            }

            currentMatRootIndex += K;
        }

    }

    SparseMatrix<complex<double>> H(pow(2, N), pow(2, N));
    H.setFromTriplets(elementList.begin(), elementList.end());
    H.diagonal().array() += N*0.75;
    H.makeCompressed();

    return H;
}

SparseMatrix<complex<double>> spinOpS2_momentum_sparse_m(int N, int n_up) {
    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    // loop over all magnetizations m
    //int n_up = N/2;
    // Find states compatible with m and store them in list.
    vector<int> s_vector_m = getStates_m(N, n_up);
    int M = s_vector_m.size();

    typedef Eigen::Triplet<complex<double>> matFil;

    // "global" index corresponding to the (1,1)-Element of the current block
    int currentMatRootIndex = 0;

    // init list of elements to be filled
    vector<matFil> elementList;

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
                for (int j = 0; j < i; j++) {
                    if (getBit(a, i) == getBit(a, j)) {
                        elementList.emplace_back(matFil(l + currentMatRootIndex, l + currentMatRootIndex, 0.5));
                    } else {
                        elementList.emplace_back(matFil(l + currentMatRootIndex, l + currentMatRootIndex, -0.5));

                        int b = a;
                        flipBit(b, i);
                        flipBit(b, j);
                        vector<int> r_L = representative(b, N);
                        int f = findState(s_vector_k, r_L[0]);

                        if (f >= 0) {
                            complex<double> offDiagEl = 2 * 0.5 * sqrt((double) R_vector[l] / (double) R_vector[f])
                                                        * std::exp(complex<double>(0, 4.0 * M_PI
                                                                                      * (double) k * (double) r_L[1] / (double) N));
                            elementList.emplace_back(matFil(l + currentMatRootIndex, f + currentMatRootIndex, offDiagEl));
                        }
                    }
                }
            }
        }
        currentMatRootIndex += K;
    }

    long siz = fact(N) / (fact(N-n_up) * fact(n_up));
    SparseMatrix<complex<double>> H(siz,siz);
    H.setFromTriplets(elementList.begin(), elementList.end());
    H.diagonal().array() += N*0.75;
    H.makeCompressed();

    return H;
}

vector<SparseMatrix<complex<double>>> momentumHamiltonian_sparse_blocks(double J_ratio, int N, int m_lower, int m_upper) {
    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    // init list of blocks
    list<SparseMatrix<complex<double>>> H_subspace_list;

    // loop over all magnetizations m
    for (int m_setter = m_lower; m_setter <= m_upper; m_setter++) {
        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, m_setter);
        int M = s_vector_m.size();

        int currentMatRootIndex = 0;
        // init (m,k)-block
        list<Eigen::Triplet<complex<double>>> elementList;

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
                    int j = (i + 1) % N;
                    if (getBit(a, i) == getBit(a, j)) {
                        elementList.emplace_back(l + currentMatRootIndex, l + currentMatRootIndex, 0.25);
                    } else {
                        elementList.emplace_back( + currentMatRootIndex, l + currentMatRootIndex, -0.25);

                        int b = a;
                        flipBit(b, i);
                        flipBit(b, j);
                        vector<int> r_L = representative(b, N);
                        int f = findState(s_vector_k, r_L[0]);

                        if (f >= 0) {
                            complex<double> offDiagEl = 0.5 * sqrt((double) R_vector[l] / (double) R_vector[f])
                                                        * std::exp(
                                    complex<double>(0, 4.0 * M_PI * (double) k * (double) r_L[1] / (double) N));
                            elementList.emplace_back(l + currentMatRootIndex, f + currentMatRootIndex, offDiagEl);
                        }
                    }
                    int n = (i + 2) % N;
                    if (getBit(a, i) == getBit(a, n) && i < N) {
                        elementList.emplace_back(l + currentMatRootIndex, l + currentMatRootIndex, J_ratio * 0.25);
                    } else if (i < N) {
                        elementList.emplace_back(l + currentMatRootIndex, l + currentMatRootIndex, -J_ratio * 0.25);

                        int b = a;
                        flipBit(b, i);
                        flipBit(b, n);
                        vector<int> r_L = representative(b, N);
                        int f = findState(s_vector_k, r_L[0]);

                        if (f >= 0) {
                            complex<double> offDiagEl = 0.5 * sqrt((double) R_vector[l] / (double) R_vector[f])
                                                        * std::exp(
                                    complex<double>(0, 4.0 * M_PI * (double) k * (double) r_L[1] / (double) N));
                            elementList.emplace_back(l + currentMatRootIndex, f + currentMatRootIndex, J_ratio * offDiagEl);
                        }
                    }
                }
            }
            currentMatRootIndex += K;
        }
        SparseMatrix<complex<double>> H_block(M, M);
        H_block.setFromTriplets(elementList.begin(), elementList.end());
        H_subspace_list.push_back(H_block);
    }
    return vector<SparseMatrix<complex<double>>>(H_subspace_list.begin(), H_subspace_list.end());
}
