#include "magnetizationHamiltonian_sparse.h"

using std::vector;
using std::list;
using Eigen::SparseMatrix;

vector<SparseMatrix<double>> magnetizationHamiltonian_sparse(double J_ratio, int N) {
    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    typedef Eigen::Triplet<double> matFil;

    // init empty block list
    vector<SparseMatrix<double>> H_subspace_list;

    // loop over all magnetizations
    for (int m_setter = 0; m_setter <= N; m_setter++) {
        // Calculate magnetization and number of "up"-states for given magnetization.
        double m = -N/2.0 + m_setter;
        int n_up = round( m + N/2.0);

        // init list of elements to be filled
        list<matFil> elementList;

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, n_up);
        int M = s_vector_m.size();

        // Generate Block with correct size and fill elements.

        for (int k = 0; k < M; k++) {
            int a = s_vector_m[k];
            for (int i = 0; i <N; i++) {
                int j = (i+1) % N;
                if (getBit(a, i) == getBit(a, j)) {
                    elementList.emplace_back(matFil(k, k, 0.25));
                } else {
                    elementList.emplace_back(matFil(k, k, -0.25));

                    int b = a;
                    flipBit(b, i);
                    flipBit(b, j);
                    int l = findState(s_vector_m, b);
                    elementList.emplace_back(matFil(k, l, 0.5));
                }
            }
            for (int i = 0; i < N; i++) {
                int j = (i+2) % N;
                if (getBit(a, i) == getBit(a, j)) {
                    elementList.emplace_back(matFil(k, k, J_ratio*0.25));
                } else {
                    elementList.emplace_back(matFil(k, k, J_ratio*-0.25));

                    int b = a;
                    flipBit(b, i);
                    flipBit(b, j);
                    int l = findState(s_vector_m, b);
                    elementList.emplace_back(matFil(k, l, J_ratio*0.5));
                }
            }
        }
        SparseMatrix<double> H_block(M, M);
        H_block.setFromTriplets(elementList.begin(), elementList.end());
        H_block.makeCompressed();
        H_subspace_list.emplace_back(H_block);
    }
    return H_subspace_list;
}

vector<SparseMatrix<double>> spinOp2_magnetization_sparse(int N) {
    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    typedef Eigen::Triplet<double> matFil;

    // init empty block list
    vector<SparseMatrix<double>> H_subspace_list;

    // loop over all magnetizations
    for (int m_setter = 0; m_setter <= N; m_setter++) {
        // Calculate magnetization and number of "up"-states for given magnetization.
        double m = -N/2.0 + m_setter;
        int n_up = round( m + N/2.0);

        // init list of elements to be filled
        list<matFil> elementList;

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, n_up);
        int M = s_vector_m.size();

        // Generate Block with correct size and fill elements.

        for (int k = 0; k < M; k++) {
            int a = s_vector_m[k];
            for (int i = 0; i <N; i++) {
                for (int j = 0; j < i; j++) {
                    if (getBit(a, i) == getBit(a, j)) {
                        elementList.emplace_back(matFil(k, k, 0.5));
                    } else {
                        elementList.emplace_back(matFil(k, k, -0.5));

                        int b = a;
                        flipBit(b, i);
                        flipBit(b, j);
                        int l = findState(s_vector_m, b);
                        elementList.emplace_back(matFil(k, l, 1));
                    }
                }
            }
        }
        for (int i = 0; i < M; i++) {
            elementList.emplace_back(i, i, N*0.75);
        }
        SparseMatrix<double> H_block(M, M);
        H_block.setFromTriplets(elementList.begin(), elementList.end());
        H_block.makeCompressed();
        H_subspace_list.emplace_back(H_block);
    }
    return H_subspace_list;
}