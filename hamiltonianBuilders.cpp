//
// This cpp file contains functions used to generate the Hamiltonian using different approaches, as well as the spin
// operator S². Also included: Functions to get eigenvalues and generate full matrices from blocks.
//

#include "hamiltonianBuilders.h"

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

// Get vector containing all states corresponding to a given m.
vector<int> getStates_m(int N, int n_up) {
    vector<int> s_vector_m;
    for (int s = 0; s < pow(2, N); s++) {
        if (bitSum(s) == n_up) {s_vector_m.push_back(s);}
    }
    return s_vector_m;
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

list<list<list<MatrixXd>>> parityHamiltonian(double J_ratio, int N) {
    // N must be even and > 6 or this no longer describes the correct system.
    //if (N < 6 || N%2 == 1) {
    //    throw std::invalid_argument("N must be larger than 6 and even.");
    //}

    // init list of block lists
    list<list<list<MatrixXd>>> H_subspace_list;
    int g = 0;

    // loop over all magnetizations m
    for (int n_up = 0; n_up <= N; n_up++) {

        int y = 0;
        // Calculate magnetization and number of "up"-states for given magnetization.

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, n_up);

        // init list of blocks
        list<list<MatrixXd>> H_subSubspace_list;

        // loop over all possible semi-momenta k and parity numbers p = +-1
        for (int k = 0; k <= trunc(N/4); k++) {
            list<MatrixXd> H_subSubSubspace_list;
            for (int p : {-1, 1}) {
                vector<int> s_vector_k, R_vector, m_vector;
                for (int s: s_vector_m) {
                    for (int sigma : {-1, 1}) {
                        if ((k == 0 || k == trunc(N/4) ) && sigma == -1) continue; //
                        vector<int> R_m = checkState_parity(s, k, N);
                        if (R_m[1] != -1) {
                            complex<double> v = sigma * (double) p * std::cos(
                                    (double) k * (double) R_m[1] * 4.0 * M_PI
                                    / (double) N);
                            if (std::abs( 1.0 + v ) < epsilon ) R_m[0] = -1;
                            if (sigma == -1 && std::abs( 1.0 - v ) > epsilon ) R_m[0] = -1;
                        }
                        if (R_m[0] > 0) {
                            s_vector_k.emplace_back(s);
                            R_vector.emplace_back(sigma * R_m[0]);
                            m_vector.emplace_back(R_m[1]);
                        }
                    }
                }
                int K = s_vector_k.size();

                MatrixXd H = MatrixXd::Zero(K, K);

                for (int a = 0; a < K; a++) {
                    const int s = s_vector_k[a];
                    g++;
                    y++;
                    //std::cout << s << " " << mag << " " << k << " " << p << " " << R_vector[a] << std::endl;
                    int n;
                    if (a > 0 && s_vector_k[a] == s_vector_k[a-1]) continue;
                    if (a < K-1 && s_vector_k[a] == s_vector_k[a+1]) {n = 2;
                    } else n = 1;

                    for (int u = a; u < a + n; u++) {
                        H(u, u) += E_z_parity(s, J_ratio, N);
                    }

                    for (int i = 0; i < N; i++) {
                        int s_prime = s;
                        int j = (i + 1) % N;
                        if (getBit(s_prime, i) != getBit(s_prime, j)) {
                            flipBit(s_prime, i);
                            flipBit(s_prime, j);
                            vector<int> r_l_q = representative_parity(s_prime, N);
                            int b = findState(s_vector_k, r_l_q[0]);
                            int m;
                            if (b >= 0) {
                                if (b > 0 && s_vector_k[b] == s_vector_k[b - 1]) {
                                    m = 2;
                                    b += -1;
                                } else if (b < K - 1 && s_vector_k[b] == s_vector_k[b + 1]) {
                                    m = 2;
                                } else m = 1;
                                for (int i_mat = a; i_mat < a + n; i_mat++) {
                                    for (int j_mat = b; j_mat < b + m; j_mat++) {
                                        double val = h_Element_parity(i_mat, j_mat, r_l_q[1], r_l_q[2], k, p, N,
                                                                      s_vector_k, R_vector, m_vector);
                                        H(i_mat, j_mat) += val;
                                    }
                                }
                            }
                        }
                    }

                    /*for (int i = 0; i < N; i++) {
                        int s_prime = s;
                        int j = (i + 2) % N;
                        if (getBit(s_prime, i) != getBit(s_prime, j)) {
                            flipBit(s_prime, i);
                            flipBit(s_prime, j);
                            vector<int> r_l_q = representative_parity(s_prime, N);
                            int b = findState(s_vector_k, r_l_q[0]);
                            int m;
                            if (b >= 0) {
                                if (b > 0 && s_vector_k[b] == s_vector_k[b - 1]) {
                                    m = 2;
                                    b += -1;
                                } else if (b < K - 1 && s_vector_k[b] == s_vector_k[b + 1]) {
                                    m = 2;
                                } else m = 1;
                                for (int j_mat = b; j_mat < b + m; j_mat++) {
                                    for (int i_mat = a; i_mat < a + n; i_mat++) {
                                        double val = h_Element_parity(i_mat, j_mat, r_l_q[1], r_l_q[2], k, p, N,
                                                                      s_vector_k, R_vector, m_vector);
                                        H(i_mat, j_mat) += J_ratio * val;

                                    }
                                }
                            }
                        }
                    }*/
                }
                H_subSubSubspace_list.emplace_back(H);
            }
            H_subSubspace_list.emplace_back(H_subSubSubspace_list);
        }
        H_subspace_list.emplace_back(H_subSubspace_list);
        //std::cout << y << "\n";
    }
    //std::cout << g << "\n";
    return H_subspace_list;
}

double h_Element_parity(int a, int b, double l, double q, double k, double p, double N,
                        const vector<int> & s_vec, const vector<int> & R_vec, const vector<int> & m_vec) {
    double sigma_a = (double) R_vec[a]/abs(R_vec[a]);
    double sigma_b = (double) R_vec[b]/abs(R_vec[b]);
    double Na = N_a_sigma(g_k(k, N), N, R_vec[a], p, k, m_vec[a]);
    double Nb = N_a_sigma(g_k(k, N), N, R_vec[b], p, k, m_vec[b]);
    double k_actual = (double) k*4.0*M_PI/ (double) N;
    double ret = 0;
    if (sigma_a == sigma_b) {
        if (m_vec[b] == -1) {
            ret = 0.5 * pow(sigma_a * p, q) * sqrt(Nb/Na) * cos(k_actual*l);
        } else {
            ret = 0.5 * pow(sigma_a * p, q) * sqrt(Nb/Na) * (cos(k_actual*l) + sigma_a * p * cos(k_actual*(l-m_vec[b])))
                / (1.0 + sigma_a * p * cos(k_actual * m_vec[b]));
        }
    } else {
        if (m_vec[b] == -1) {
            ret = 0.5 * pow(sigma_a * p, q) * sqrt(Nb/Na) * -sigma_a * sin(k_actual*l);
        } else {
            ret = 0.5 * pow(sigma_a * p, q) * sqrt(Nb/Na) * (-sigma_a * sin(k_actual*l) + p * sin(k_actual*(l-m_vec[b])))
                / (1.0 - sigma_a * p * cos(k_actual * m_vec[b]));
        }
    }
    return ret;
}

double g_k(double k, double N) {
    if (k < epsilon || (k - trunc(N/4)) < epsilon) {
        return 2;
    } else {
        return 1;
    }
}

double N_a_sigma(double g, double N, double sigmaR, double p, double k, double m) {
    if (m == -1) {
        return N * N * g /abs(sigmaR);
    } else {
        return N * N * g / abs(sigmaR) * (1.0 + sigmaR / abs(sigmaR)
                * p * cos(k*m*4*M_PI/N));
    }
}

double E_z_parity(const int s, const double J_ratio, const int N) {
    double E_z = 0;
    for (int i = 0; i < N; i++) {
        int j = (i + 1) % N;
        if (getBit(s, i) == getBit(s, j)) {
            E_z += 0.25;
        } else {
            E_z += -0.25;
        }
        j = (i + 2) % N;
        if (getBit(s, i) == getBit(s, j)) {
            E_z += J_ratio * 0.25;
        } else {
            E_z += J_ratio * -0.25;
        }
    }
    return E_z;
}

// Diagonalization methods for a number of cases.
vector<double> getEnergiesFromBlocks(const list<list<list<MatrixXd>>> & H_list, int N) {
    vector<double> energies;
    for (const list<list<MatrixXd>> & H_sublist : H_list) {
        for (const list<MatrixXd> & H_subSubList : H_sublist) {
            for (const MatrixXd & mat : H_subSubList) {
                Eigen::SelfAdjointEigenSolver<MatrixXd> sol;
                if (mat.cols() == 0) {
                    continue;
                }
                sol.compute(mat);
                Eigen::VectorXcd blockEnergies = sol.eigenvalues();
                std::for_each(blockEnergies.begin(), blockEnergies.end(),
                              [&](complex<double> & d){energies.emplace_back(d.real());});
            }
        }
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<MatrixXcd> & H_list, int N) {
    vector<double> energies(pow(2, N), 0);
    int j = 0;
    for (const MatrixXcd & mat : H_list) {
        Eigen::SelfAdjointEigenSolver<MatrixXcd> sol;
        if (mat.cols() == 0) {
            continue;
        }
        sol.compute(mat);
        Eigen::VectorXcd blockEnergies = sol.eigenvalues();
        std::for_each(blockEnergies.begin(), blockEnergies.end(),
                      [&](complex<double> & d){energies[j] = d.real(); j++;});
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<MatrixXcd> & H_list) {
    vector<double> energies;
    for (const MatrixXcd & mat : H_list) {
        Eigen::SelfAdjointEigenSolver<MatrixXcd> sol;
        if (mat.cols() == 0) {
            continue;
        }
        sol.compute(mat);
        Eigen::VectorXcd blockEnergies = sol.eigenvalues();
        std::for_each(blockEnergies.begin(), blockEnergies.end(),
                      [&](complex<double> & d){energies.emplace_back( d.real() );});
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<MatrixXd> & H_list, int N) {
    vector<double> energies(pow(2, N), 0);
    int j = 0;
    for (const MatrixXd & mat : H_list) {
        Eigen::SelfAdjointEigenSolver<MatrixXd> sol;
        if (mat.cols() == 0) {
            continue;
        }
        sol.compute(mat);
        Eigen::VectorXcd blockEnergies = sol.eigenvalues();
        std::for_each(blockEnergies.begin(), blockEnergies.end(),
                      [&](complex<double> & d){energies[j] = d.real(); j++;});
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<MatrixXd> & H_list) {
    vector<double> energies;
    for (const MatrixXd & mat : H_list) {
        Eigen::SelfAdjointEigenSolver<MatrixXd> sol;
        if (mat.cols() == 0) {
            continue;
        }
        sol.compute(mat);
        Eigen::VectorXcd blockEnergies = sol.eigenvalues();
        std::for_each(blockEnergies.begin(), blockEnergies.end(),
                      [&](complex<double> & d){energies.emplace_back( d.real() );});
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<list<MatrixXcd>> & H_list, int N) {
    vector<double> energies(pow(2, N), 0);
    int j = 0;
    for (const list<MatrixXcd> & subList : H_list) {
        for (const MatrixXcd & mat : subList) {
            Eigen::SelfAdjointEigenSolver<MatrixXcd> sol;
            if (mat.cols() == 0) {
                continue;
            }
            sol.compute(mat);
            Eigen::VectorXcd blockEnergies = sol.eigenvalues();
            std::for_each(blockEnergies.begin(), blockEnergies.end(),
                          [&](complex<double> &d) {
                              energies[j] = d.real();
                              j++; });
        }
    }
    std::sort(energies.begin(), energies.end());
    return energies;
}



vector<vector<vector<double>>> getEnergiesFromBlocksByK(const list<list<MatrixXcd>> & H_list) {
    vector<vector<vector<double>>> energies;

    for (const list<MatrixXcd> & subList : H_list) {
        vector<vector<double>> m_ergs;
        for (const MatrixXcd & mat : subList) {
            Eigen::SelfAdjointEigenSolver<MatrixXcd> sol;
            if (mat.cols() == 0) {
                continue;
            }
            sol.compute(mat);
            Eigen::VectorXd blockEnergies = sol.eigenvalues().real();
            vector<double> k_ergs(blockEnergies.begin(), blockEnergies.end());
            m_ergs.emplace_back(k_ergs);
        }
        energies.emplace_back(m_ergs);
    }
    return energies;
}


// Generate full-sized matrix from blocks.
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
    for (const list<MatrixXcd> & l : matrix_doubleList) {
        int size = 0;
        for (const MatrixXcd & m : l) {
            size += m.rows();
        }
        bdmlist.emplace_back( blkdiag(l, size) );
    }
    return bdmlist;
}

// Save vector containing eigenvalues to file.
void saveEnergies(const vector<double> & ergs, const std::string & path) {
    std::ofstream ergFile;
    ergFile.open(path);
    for (double d : ergs) {
        ergFile << d << "\n";
    }
    ergFile.close();
}