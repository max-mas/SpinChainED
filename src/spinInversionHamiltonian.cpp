//
// Created by MaxM on 16/05/2022.
//

#include "spinInversionHamiltonian.h"

using std::list;
using std::vector;

using Eigen::MatrixXd;

list<list<list<MatrixXd>>> spinInversionHamiltonian(double J_ratio, int N, int n_up_min, int n_up_max) {
    // N must be a multiple of 4 and >=8.
    if (N < 8 || N%4 != 0) {
        throw std::invalid_argument("N must be larger than 8 and a multiple of 4.");
    }

    // init list of block lists
    list<list<list<MatrixXd>>> H_subspace_list;

    // loop over all magnetizations m
    for (int n_up = n_up_min; n_up <= n_up_max; n_up++) {

        // Calculate magnetization and number of "up"-states for given magnetization.

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, n_up);

        // init list of blocks
        list<list<MatrixXd>> H_subSubspace_list;

        if (n_up == N/2) {
            for (int k = 0; k <= trunc(N/4); k++) {
                list<MatrixXd> H_subSubSubspace_list;
                for (int p : {-1, 1}) {
                    for (int z : {-1, 1}) {

                        vector<int> s_vector_k, R_vector, m_vector, n_vector, c_vector;

                        getStates_k_p_z(N, s_vector_m, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector);

                        int K = s_vector_k.size();
                        MatrixXd H = MatrixXd::Zero(K, K);

                        for (int a = 0; a < K; a++) {
                            const int s = s_vector_k[a];
                            int n;
                            if (a > 0 && s_vector_k[a] == s_vector_k[a-1]) continue;
                            if (a < K-1 && s_vector_k[a] == s_vector_k[a+1]) {n = 2;
                            } else n = 1;

                            for (int u = a; u < a + n; u++) {
                                H(u, u) += E_z_parity(s, J_ratio, N);
                            }

                            setHElementsForState_inversion(N, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector, K, H, a, s, n, 1, 1);
                            setHElementsForState_inversion(N, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector, K, H, a, s, n, 2, J_ratio);
                        }
                        H_subSubSubspace_list.emplace_back(H);
                    }
                }
                H_subSubspace_list.emplace_back(H_subSubSubspace_list);
            }

            H_subspace_list.emplace_back(H_subSubspace_list);
            continue;
        }

        for (int k = 0; k <= trunc(N/4); k++) {
            list<MatrixXd> H_subSubSubspace_list;
            for (int p : {-1, 1}) {

                vector<int> s_vector_k, R_vector, m_vector;
                getStates_k_p(N, s_vector_m, k, p, s_vector_k, R_vector, m_vector);
                int K = s_vector_k.size();

                MatrixXd H = MatrixXd::Zero(K, K);

                for (int a = 0; a < K; a++) {
                    const int s = s_vector_k[a];
                    int n;
                    if (a > 0 && s_vector_k[a] == s_vector_k[a-1]) continue;
                    if (a < K-1 && s_vector_k[a] == s_vector_k[a+1]) {n = 2;
                    } else n = 1;

                    for (int u = a; u < a + n; u++) {
                        H(u, u) += E_z_parity(s, J_ratio, N);
                    }

                    setHElementsForState_parity(N, k, p, s_vector_k, R_vector, m_vector, K, H, a, s, n, 1, 1);
                    setHElementsForState_parity(N, k, p, s_vector_k, R_vector, m_vector, K, H, a, s, n, 2, J_ratio);

                }
                H_subSubSubspace_list.emplace_back(H);
            }
            H_subSubspace_list.emplace_back(H_subSubSubspace_list);
        }
        H_subspace_list.emplace_back(H_subSubspace_list);

    }
    return H_subspace_list;
}

vector<double> getEnergies_memorySaving_threaded_inversion(double J_ratio, int N) {
    // N must be a multiple of 4 and >=8.
    if (N < 8 || N%4 != 0) {
        throw std::invalid_argument("N must be larger than 8 and a multiple of 4.");
    }

    // init list of block lists
    vector<double> ergs;

    // loop over all magnetizations m
#pragma omp parallel for default(none) shared(ergs, J_ratio, N)
    for (int n_up = 0; n_up <= N; n_up++) {

        // Calculate magnetization and number of "up"-states for given magnetization.

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, n_up);

        if (n_up == N/2) {
            for (int k = 0; k <= trunc(N/4); k++) {
                for (int p = -1; p <= 1; p+= 2) {
                    for (int z = -1; z <= 1; z+= 2) {

                        vector<int> s_vector_k, R_vector, m_vector, n_vector, c_vector;

                        getStates_k_p_z(N, s_vector_m, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector);

                        int K = s_vector_k.size();
                        MatrixXd H = MatrixXd::Zero(K, K);

                        for (int a = 0; a < K; a++) {
                            const int s = s_vector_k[a];
                            int n;
                            if (a > 0 && s_vector_k[a] == s_vector_k[a-1]) continue;
                            if (a < K-1 && s_vector_k[a] == s_vector_k[a+1]) {n = 2;
                            } else n = 1;

                            for (int u = a; u < a + n; u++) {
                                H(u, u) += E_z_parity(s, J_ratio, N);
                            }

                            setHElementsForState_inversion(N, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector, K, H, a, s, n, 1, 1);
                            setHElementsForState_inversion(N, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector, K, H, a, s, n, 2, J_ratio);
                        }

                        Eigen::SelfAdjointEigenSolver<MatrixXd> sol;
                        if (H.cols() == 0) {
                            continue;
                        }
                        sol.compute(H);
                        Eigen::VectorXd blockErgs = sol.eigenvalues().real();
                        vector<double> blockErgs_stl(blockErgs.begin(), blockErgs.end());
                        writeThreadSafe(ergs, blockErgs_stl);
                    }
                }
            }
            continue;
        }

        for (int k = 0; k <= trunc(N/4); k++) {
            for (int p = -1; p <= 1; p++) {

                vector<int> s_vector_k, R_vector, m_vector;
                getStates_k_p(N, s_vector_m, k, p, s_vector_k, R_vector, m_vector);
                int K = s_vector_k.size();

                MatrixXd H = MatrixXd::Zero(K, K);

                for (int a = 0; a < K; a++) {
                    const int s = s_vector_k[a];
                    int n;
                    if (a > 0 && s_vector_k[a] == s_vector_k[a-1]) continue;
                    if (a < K-1 && s_vector_k[a] == s_vector_k[a+1]) {n = 2;
                    } else n = 1;

                    for (int u = a; u < a + n; u++) {
                        H(u, u) += E_z_parity(s, J_ratio, N);
                    }

                    setHElementsForState_parity(N, k, p, s_vector_k, R_vector, m_vector, K, H, a, s, n, 1, 1);
                    setHElementsForState_parity(N, k, p, s_vector_k, R_vector, m_vector, K, H, a, s, n, 2, J_ratio);

                }
                Eigen::SelfAdjointEigenSolver<MatrixXd> sol;
                if (H.cols() == 0) {
                    continue;
                }
                sol.compute(H);
                Eigen::VectorXd blockErgs = sol.eigenvalues().real();
                vector<double> blockErgs_stl(blockErgs.begin(), blockErgs.end());
                writeThreadSafe(ergs, blockErgs_stl);
            }
        }
    }
    std::sort(ergs.begin(), ergs.end());
    return ergs;
}

list<list<MatrixXd>> spinOpS2_spinInv_m0(int N) {
    vector<int> s_vector_m = getStates_m(N, N/2);

    // init list of blocks
    list<list<MatrixXd>> S2_subSubspace_list;

    for (int k = 0; k <= trunc(N/4); k++) {
        list<MatrixXd> S2_subSubSubspace_list;
        for (int p : {-1, 1}) {
            for (int z : {-1, 1}) {

                vector<int> s_vector_k, R_vector, m_vector, n_vector, c_vector;

                getStates_k_p_z(N, s_vector_m, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector);

                int K = s_vector_k.size();
                MatrixXd S2 = N*0.75*MatrixXd::Identity(K, K);

                for (int a = 0; a < K; a++) {
                    const int s = s_vector_k[a];
                    int n;
                    if (a > 0 && s_vector_k[a] == s_vector_k[a-1]) continue;
                    if (a < K-1 && s_vector_k[a] == s_vector_k[a+1]) {n = 2;
                    } else n = 1;

                    for (int u = a; u < a + n; u++) {
                        double E_z = 0;
                        for (int i = 0; i < N; i++) {
                            for (int j = 0; j < i; j++) {
                                if (getBit(s, i) == getBit(s, j)) {
                                    E_z += 0.25;
                                } else {
                                    E_z += -0.25;
                                }
                            }

                        }
                        S2(u, u) += 2*E_z_parity(s, 1, N);
                    }
                    for (int i = 0; i < N; i++) {
                        int s_prime = s;
                        for (int j = 0; j < i; j++) {
                            if (getBit(s_prime, i) != getBit(s_prime, j)) {
                                flipBit(s_prime, i);
                                flipBit(s_prime, j);
                                vector<int> r_l_q_g = representative_inversion(s_prime, N);
                                int b = findState(s_vector_k, r_l_q_g[0]);
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
                                            double val = 2*hElement_inversion(i_mat, j_mat, r_l_q_g[1], r_l_q_g[2], r_l_q_g[3],
                                                                              k, p, z, N, R_vector, m_vector, n_vector, c_vector);
                                            S2(i_mat, j_mat) += val;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                S2_subSubSubspace_list.emplace_back(S2);
            }
        }
        S2_subSubspace_list.emplace_back(S2_subSubSubspace_list);
    }
    return S2_subSubspace_list;
}


void
getStates_k_p_z(int N, vector<int> &s_vector_m, int k, int p, int z, vector<int> &s_vector_k, vector<int> &R_vector,
                vector<int> &m_vector, vector<int> &n_vector, vector<int> &c_vector) {
    for (int s : s_vector_m) {
        for (int sigma : {-1, 1}) {
            if ((k == 0 || k == trunc(N/4) ) && sigma == -1) continue;

            vector<int> R_mp_mz_mpz = checkState_inversion(s, k, N);
            int m = -1;
            int n = -1;
            int c = getStateClass(R_mp_mz_mpz[1], R_mp_mz_mpz[2], R_mp_mz_mpz[3]);

            if (c == 2) {
                m = R_mp_mz_mpz[1];
                double Na = N_a_sigma_inv(N, sigma * R_mp_mz_mpz[0], k, p, z, m, n, c);
                double Na_min = N_a_sigma_inv(N, -sigma * R_mp_mz_mpz[0], k, p, z, m, n, c);
                if (abs(Na) < epsilon) R_mp_mz_mpz[0] = -1;
                if (sigma == -1 && Na_min > epsilon) R_mp_mz_mpz[0] = -1;
            } else if (c == 3) {
                m = R_mp_mz_mpz[2];
                double Na = N_a_sigma_inv(N, sigma * R_mp_mz_mpz[0], k, p, z, m, n, c);
                if (abs(Na) < epsilon) R_mp_mz_mpz[0] = -1;
            } else if (c == 4) {
                m = R_mp_mz_mpz[3];
                double Na = N_a_sigma_inv(N, sigma * R_mp_mz_mpz[0], k, p, z, m, n, c);
                double Na_min = N_a_sigma_inv(N, -sigma * R_mp_mz_mpz[0], k, p, z, m, n, c);
                if (abs(Na) < epsilon) R_mp_mz_mpz[0] = -1;
                if (sigma == -1 && Na_min > epsilon) R_mp_mz_mpz[0] = -1;
            } else if (c == 5) {
                m = R_mp_mz_mpz[1];
                n = R_mp_mz_mpz[2];
                double Na = N_a_sigma_inv(N, sigma * R_mp_mz_mpz[0], k, p, z, m, n, c);
                double Na_min = N_a_sigma_inv(N, -sigma * R_mp_mz_mpz[0], k, p, z, m, n, c);
                if (abs(Na) < epsilon) R_mp_mz_mpz[0] = -1;
                if (sigma == -1 && Na_min > epsilon) R_mp_mz_mpz[0] = -1;
            }

            if (R_mp_mz_mpz[0] > 0) {
                s_vector_k.emplace_back(s);
                R_vector.emplace_back(sigma * R_mp_mz_mpz[0]);
                m_vector.emplace_back(m);
                n_vector.emplace_back(n);
                c_vector.emplace_back(c);
            }
        }
    }
}

void setHElementsForState_inversion(const int N, const int k, const int p, const int z, const vector<int> & s_vec,
                                    const vector<int> & R_vec, const vector<int> & m_vec, const vector<int> & n_vec,
                                    const vector<int> & c_vec, const int K, MatrixXd & H, const int a, const int s,
                                    const int n, const int neighbour, const double J) {
    for (int i = 0; i < N; i++) {
        int s_prime = s;
        int j = (i + neighbour) % N;
        if (getBit(s_prime, i) != getBit(s_prime, j)) {
            flipBit(s_prime, i);
            flipBit(s_prime, j);
            vector<int> r_l_q_g = representative_inversion(s_prime, N);
            int b = findState(s_vec, r_l_q_g[0]);
            int m;
            if (b >= 0) {
                if (b > 0 && s_vec[b] == s_vec[b - 1]) {
                    m = 2;
                    b += -1;
                } else if (b < K - 1 && s_vec[b] == s_vec[b + 1]) {
                    m = 2;
                } else m = 1;
                for (int i_mat = a; i_mat < a + n; i_mat++) {
                    for (int j_mat = b; j_mat < b + m; j_mat++) {
                        double val = hElement_inversion(i_mat, j_mat, r_l_q_g[1], r_l_q_g[2], r_l_q_g[3],
                                                            k, p, z, N, R_vec, m_vec, n_vec, c_vec);
                        H(i_mat, j_mat) += J * val;
                    }
                }
            }
        }
    }
}

double hElement_inversion(const int a, const int b, const double l, const double q, const double g, const double k, const double p, const double z, const double N,
                          const vector<int> & R_vec, const vector<int> & m_vec, const vector<int> & n_vec,
                          const vector<int> & c_vec) {
    double sigma_a = (double) R_vec[a] / (double) abs(R_vec[a]);
    double sigma_b = (double) R_vec[b] / (double) abs(R_vec[b]);
    double Na = N_a_sigma_inv(N, R_vec[a], k, p, z, m_vec[a], n_vec[a], c_vec[a]);
    double Nb = N_a_sigma_inv(N, R_vec[b], k, p, z, m_vec[b], n_vec[b], c_vec[b]);
    double k_actual = (double) k * 4.0 * M_PI / (double) N;
    double ret = 0.5 * pow(sigma_a * p, q) * pow(z, g) * sqrt(Nb/Na);

    if (sigma_a == sigma_b) {
        if (c_vec[b] == 1 || c_vec[b] == 3) {
            ret *= cos(k_actual * l);
        } else if (c_vec[b] == 2 || c_vec[b] == 5) {
            ret *= (cos(k_actual * l) + sigma_a * p * cos(k_actual*(l - m_vec[b]))) / (1.0 + sigma_a * p * cos(k_actual * m_vec[b]));
        } else if (c_vec[b] == 4) {
            ret *= (cos(k_actual * l) + sigma_a * p * z * cos(k_actual*(l - m_vec[b]))) / (1.0 + sigma_a * p * z * cos(k_actual * m_vec[b]));
        }
    } else {
        if (c_vec[b] == 1 || c_vec[b] == 3) {
            ret *= -sin(k_actual * l);
        } else if (c_vec[b] == 2 || c_vec[b] == 5) {
            ret *= (-sin(k_actual * l) + p * sin(k_actual*(l - m_vec[b]))) / (1.0 - sigma_a * p * cos(k_actual * m_vec[b]));
        } else if (c_vec[b] == 4) {
            ret *= (-sin(k_actual * l) + p * z * sin(k_actual*(l - m_vec[b]))) / (1.0 - sigma_a * p * z * cos(k_actual * m_vec[b]));
        }
    }
    return ret;
}

double N_a_sigma_inv(const double N, const double sigmaR, const double k, const double p, const double z, const double m, const double n, const int c) {
    double k_actual = k * 4.0 * M_PI / N;
    double ret = 2.0 * pow(N, 2) / (abs(sigmaR) * g_k(k, N));
    double sigma = sigmaR/abs(sigmaR);
    if (c == 1) {
        ret *= 1;
    } else if (c == 2) {
        ret *= 1.0 + sigma * p * cos(k_actual * m);
    } else if (c == 3) {
        ret *= 1.0 + z * cos(k_actual * m);
    } else if (c == 4) {
        ret *= 1.0 + sigma * p * z * cos(k_actual * m);
    } else if (c == 5) {
        ret *= (1.0 + sigma * p * cos(k_actual * m)) * (1.0 + z * cos(k_actual * n));
    }
    return ret;
}

int getStateClass(const int m_p, const int m_z, const int m_pz) {
    if (m_p == -1 && m_z == -1 && m_pz == -1) {
        return 1;
    } else if (m_p != -1 && m_z == -1 && m_pz == -1) {
        return 2;
    } else if (m_p == -1 && m_z != -1 && m_pz == -1) {
        return 3;
    } else if (m_p == -1 && m_z == -1 && m_pz != -1) {
        return 4;
    } else {
        return 5;
    }
}


vector<int> checkState_inversion(const int s, const int k, const int N) {
    int R = -1;
    int m_p = -1;
    int m_z = -1;
    int m_pz = -1;

    int t = s;
    for (int i = 1; i <= N/2; i++) {
        cycleBits2(t, N); //translate state
        if (t < s) {return {R, m_p, m_z, m_pz};}
        else if (t == s) {
            if ( k % (int) trunc(N/(2*i)) ) {return {R, m_p, m_z, m_pz};} // check compatibility with k
            R = i;
            break;
        }
    }

    t = s;
    reflectBits(t, N);
    for (int i = 0; i < R; i++) {
        if (t < s) {
            R = -1;
            return {R, m_p, m_z, m_pz};
        } else if (t == s) {
            m_p = i;
            break;
        }
        cycleBits2(t, N);
    }

    t = s;
    invertBits(t, N);
    for (int i = 0; i < R; i++) {
        if (t < s) {
            R = -1;
            return {R, m_p, m_z, m_pz};
        } else if (t == s) {
            m_z = i;
            break;
        }
        cycleBits2(t, N);
    }

    t = s;
    invertBits(t, N);
    reflectBits(t, N);
    for (int i = 0; i < R; i++) {
        if (t < s) {
            R = -1;
            return {R, m_p, m_z, m_pz};
        } else if (t == s) {
            m_pz = i;
            break;
        }
        cycleBits2(t, N);
    }

    return {R, m_p, m_z, m_pz};
}

vector<int> representative_inversion(const int s, const int N) {
    int l = 0;
    int q = 0;
    int g = 0;

    int r = s;
    int t = s;
    for (int i = 1; i <= N/2; i++) {
        cycleBits2(t, N);
        if (t < r) {
            r = t;
            l = i;
        }
    }

    t = s;
    reflectBits(t, N);
    for (int i = 0; i <= N/2; i++) {
        if (t < r) {
            r = t;
            l = i;
            q = 1;
        }
        cycleBits2(t, N);
    }

    t = s;
    invertBits(t, N);
    for (int i = 0; i <= N/2; i++) {
        if (t < r) {
            r = t;
            l = i;
            q = 0;
            g = 1;
        }
        cycleBits2(t, N);
    }

    t = s;
    invertBits(t, N);
    reflectBits(t, N);
    for (int i = 0; i <= N/2; i++) {
        if (t < r) {
            r = t;
            l = i;
            q = 1;
            g = 1;
        }
        cycleBits2(t, N);
    }

    return {r, l, q, g};
}

void invertBits(int & s, int N) {
    for (int i = 0; i < N; i++) {
        setBit(s, i, !getBit(s, i));
    }
}
