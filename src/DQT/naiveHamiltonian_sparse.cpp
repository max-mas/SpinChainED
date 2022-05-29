
#include "naiveHamiltonian_sparse.h"

using std::vector;

using Eigen::SparseMatrix;

SparseMatrix<double> naiveHamiltonian_sparse(double J_ratio, int N) {
    // N must be even and > 6 or this no longer describes the correct system.
    if (N < 6 || N%2 == 1) {
        throw std::invalid_argument("N must be larger than 6 and even.");
    }

    // init Matrix
    int matrixSize = pow(2, N);

    typedef Eigen::Triplet<double> matFil;

    // init list of elements to be filled
    vector<matFil> elementList;

    // set Matrix Elements for all states
    for (int a = 0; a < matrixSize; a++) {
        for (int i = 0; i < N; i++) {
            int j = (i+1) % N;
            int k = (i+2) % N;
            if (getBit(a, i) == getBit(a, j)) {
                elementList.emplace_back(a, a, 0.25);
            } else {
                elementList.emplace_back(a, a, -0.25);
                int b = a;
                flipBit(b, i);
                flipBit(b, j);
                elementList.emplace_back(a, b, 0.25);
            }
        }
        // J1 coupling
        for (int i = 0; i < N; i++) {
            int j = (i+2) % N;
            if (getBit(a, i) == getBit(a, j)) {
                elementList.emplace_back(a, a, J_ratio * 0.25);
            } else {
                elementList.emplace_back(a, a, J_ratio * -0.25);
                int b = a;
                flipBit(b, i);
                flipBit(b, j);
                elementList.emplace_back(a, b, J_ratio * 0.25);
            }
        }
    }
    SparseMatrix<double> H(pow(2, N), pow(2, N));
    H.setFromTriplets(elementList.begin(), elementList.end());
    H.makeCompressed();
    return H;
}
