//
// Created by mmaschke on 16/05/22.
//

#include "spinInversionHamiltonian.h"

using std::list;
using std::vector;

using Eigen::MatrixXd;



list<list<list<MatrixXd>>> spinInversionHamiltonian(double J_ratio, int N) {
    // N must be a multiple of 4 and >=8.
    if (N < 8 || N%4 != 0) {
        throw std::invalid_argument("N must be larger than 8 and a multiple of 4.");
    }

    // init list of block lists
    list<list<list<MatrixXd>>> H_subspace_list;

    // loop over all magnetizations m
    for (int n_up = 0; n_up <= N; n_up++) {

        // Calculate magnetization and number of "up"-states for given magnetization.

        // Find states compatible with m and store them in list.
        vector<int> s_vector_m = getStates_m(N, n_up);

        // init list of blocks
        list<list<MatrixXd>> H_subSubspace_list;

        if (n_up == N/2) {
            //// TODO
            for (int k = 0; k < trunc(N/4); k++) {
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

                            setHElementsForState_inversion(N, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector,
                                                                K, H, a, s, n, 1, 1);
                            setHElementsForState_inversion(N, k, p, z, s_vector_k, R_vector, m_vector, n_vector, c_vector,
                                                           K, H, a, s, n, 2, J_ratio);

                        }
                        H_subSubSubspace_list.emplace_back(H);
                    }
                }
                H_subSubspace_list.emplace_back(H_subSubSubspace_list);
            }
            H_subspace_list.emplace_back(H_subSubspace_list);
            continue;
        }

        // loop over all possible semi-momenta k and parity numbers p = +-1
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

void setHElementsForState_inversion(int N, int k, int p, int z, const vector<int> & s_vector_k,
                                    const vector<int> & R_vector, const vector<int> &m_vector,
                                    const vector<int> & n_vector, const vector<int> & c_vector, int K, MatrixXd & H,
                                    int a, const int s, int n, int neighbour, double J) {
    for (int i = 0; i < N; i++) {
        int s_prime = s;
        int j = (i + neighbour) % N;
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
                        double val = h_Element_inversion(i_mat, j_mat, r_l_q_g[1], r_l_q_g[2], r_l_q_g[3], k, p, z, N,
                                                      s_vector_k, R_vector, m_vector, n_vector, c_vector);
                        H(i_mat, j_mat) += J * val;
                    }
                }
            }
        }
    }
}

double h_Element_inversion(int a, int b, double l, double q, double g, double k, double p, double z, double N,
                           const vector<int> & s_vec, const vector<int> & R_vec, const vector<int> & m_vec,
                           const vector<int> & n_vec, const vector<int> & c_vec) {
    double sigma_a = (double) R_vec[a]/abs(R_vec[a]);
    double sigma_b = (double) R_vec[b]/abs(R_vec[b]);
    double k_actual = (double) k*4.0*M_PI/ (double) N;
    double Na = N_a_sigma_inversion(N, c_vec[a], sigma_a, abs(R_vec[a]), p, k_actual, m_vec[a], n_vec[a], z);
    double Nb = N_a_sigma_inversion(N, c_vec[b], sigma_b, abs(R_vec[b]), p, k_actual, m_vec[b], n_vec[b], z);
    double ret = 0;
    if (sigma_a == sigma_b) {
        if (c_vec[b] == 1 || c_vec[b] == 3) {
            ret = 0.5 * pow(sigma_a * p, q) * pow(z,g) * sqrt(Nb/Na) * cos(k_actual*l);
        } else if (c_vec[b] == 2 || c_vec[b] == 5) {
            ret = 0.5 * pow(sigma_a * p, q) * pow(z,g) * sqrt(Nb/Na) * ( cos(k_actual * l) + sigma_a * p * cos(k_actual * (l - m_vec[b])))
                    / (1 + sigma_a * p * cos(k_actual * m_vec[b]));
        } else {
            ret = 0.5 * pow(sigma_a * p, q) * pow(z,g) * sqrt(Nb/Na) * ( cos(k_actual * l) + sigma_a * p * z * cos(k_actual * (l - m_vec[b])))
                  / (1 + sigma_a * p * z * cos(k_actual * m_vec[b]));
        }
    } else {
        if (c_vec[b] == 1 || c_vec[b] == 3) {
            ret = 0.5 * pow(sigma_a * p, q) * pow(z,g) * sqrt(Nb/Na) * -sin(k_actual*l);
        } else if (c_vec[b] == 2 || c_vec[b] == 5) {
            ret = 0.5 * pow(sigma_a * p, q) * pow(z,g) * sqrt(Nb/Na) * ( -sin(k_actual * l) + p * sin(k_actual * (l - m_vec[b])))
                  / (1 - sigma_a * p * cos(k_actual * m_vec[b]));
        } else {
            ret = 0.5 * pow(sigma_a * p, q) * pow(z,g) * sqrt(Nb/Na) * ( -sin(k_actual * l) + p * z * sin(k_actual * (l - m_vec[b])))
                  / (1 - sigma_a * p * z * cos(k_actual * m_vec[b]));
        }
    }
    return ret;
}

double N_a_sigma_inversion(double N, int c, double sigma, double R, double p, double k, double z, double m, double n) {
    if (c == 1) {
        return (double) 2 * N * N / (R * g_k(k, N));
    } else if (c == 2) {
        return (double) 2 * N * N / (R * g_k(k, N)) * (1 + sigma * p * cos(k * m));
    } else if (c == 3) {
        return (double) 2 * N * N / (R * g_k(k, N)) * (1 + z * cos(k * m));
    } else if (c == 4) {
        return (double) 2 * N * N / (R * g_k(k, N)) * (1 + sigma * p * z * cos(k * m));
    } else {
        return (double) 2 * N * N / (R * g_k(k, N)) * (1 + sigma * p * cos(k * m)) * (1 + z * cos(k*n));
    }
}

void getStates_k_p_z(int N, const vector<int> &s_vector_m, int k, int p, int z, vector<int> &s_vector_k, vector<int> &R_vector,
                   vector<int> &m_vector, vector<int> & n_vector, vector<int> & c_vector) {
    for (int s: s_vector_m) {
        for (int sigma : {-1, 1}) {
            if ((k == 0 || k == trunc(N/4) ) && sigma == -1) continue; //

            vector<int> R_mp_mz_mpz = checkState_inversion(s, k, N);
            if (R_mp_mz_mpz[0] == -1) {
                continue;
            }
            int c = getStateClass(R_mp_mz_mpz[1], R_mp_mz_mpz[2], R_mp_mz_mpz[3]);
            if (c == 2 || c == 5) {
                double v = sigma * (double) p * cos(
                        (double) k * (double) R_mp_mz_mpz[1] * 4.0 * M_PI
                        / (double) N);
                if (std::abs( 1.0 + v ) < epsilon ) R_mp_mz_mpz[0] = -1;
                if (sigma == -1 && std::abs( 1.0 - v ) > epsilon ) R_mp_mz_mpz[0] = -1;
            } else if (c == 4) {
                double v = sigma * (double) p * z * cos(
                        (double) k * (double) R_mp_mz_mpz[1] * 4.0 * M_PI
                        / (double) N);
                if (std::abs( 1.0 + v ) < epsilon ) R_mp_mz_mpz[0] = -1;
                if (sigma == -1 && std::abs( 1.0 - v ) > epsilon ) R_mp_mz_mpz[0] = -1;
            }
            if (R_mp_mz_mpz[0] > 0) {
                s_vector_k.emplace_back(s);
                R_vector.emplace_back(sigma * R_mp_mz_mpz[0]);
                m_vector.emplace_back(R_mp_mz_mpz[1]);
                n_vector.emplace_back(R_mp_mz_mpz[2]);
                c_vector.emplace_back(c);
            }
        }
    }
}

vector<int> checkState_inversion(int s, int k, int N) {
    int t = s;
    int R    = -1;
    int m_p  = -1;
    int m_z  = -1;
    int m_zp = -1;
    for (int i = 1; i <= N/2; i++) {
        cycleBits2(t, N); //translate state
        if (t < s) {return {R, m_p, m_z, m_zp};}
        else if (t == s) {
            if ( k % (int) trunc(N/(2*i)) ) {return {R, m_p, m_z, m_zp};} // check compatibility with k
            R = i;
            break;
        }
    }

    t = s;
    reflectBits(t, N);
    for (int i = 0; i < R; i++) {
        if (t < s) {
            R = -1;
            return {R, m_p, m_z, m_zp};
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
            return {R, m_p, m_z, m_zp};
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
            return {R, m_p, m_z, m_zp};
        } else if (t == s) {
            m_zp = i;
            break;
        }
        cycleBits2(t, N);
    }

    return {R, m_p, m_z, m_zp};
}

vector<int> representative_inversion(const int s, const int N) {
    int r = s;
    int t = s;
    int l = 0;

    for (int i = 1; i <= N/2; i++) {
        cycleBits2(t, N);
        if (t < r) {r = t; l = i;}
    }

    //t = s;
    reflectBits(t, N);
    int q = 0;
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
    int g = 0;
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

int getStateClass(int m_p, int m_z, int m_zp) {
    if (m_p == -1 && m_z == -1 && m_zp == -1) {
        return 1;
    } else if (m_p != -1 && m_z == -1 && m_zp == -1) {
        return 2;
    } else if (m_p == -1 && m_z != -1 && m_zp == -1) {
        return 3;
    } else if (m_p == -1 && m_z == -1 && m_zp != -1) {
        return 4;
    } else {
        return 5;
    }
}


void invertBits(int & s, int N) {
        s = pow(2, N-1) - s;
}
