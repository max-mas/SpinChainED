//
// Created by MaxM on 13/05/2022.
//

#include "parityHamiltonian.h"

using std::vector;
using std::complex;
using std::list;
using std::pow;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;


list<list<list<MatrixXd>>> parityHamiltonian(double J_ratio, int N) {
    // N must be a multiple of 4 and >=8.
    if (N < 8 || N%4 != 0) {
        throw std::invalid_argument("N must be larger than 8 and a multiple of 4.");
    }

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

                    for (int i = 0; i < N; i++) {
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
                                for (int i_mat = a; i_mat < a + n; i_mat++) {
                                    for (int j_mat = b; j_mat < b + m; j_mat++) {
                                        double val = h_Element_parity(i_mat, j_mat, r_l_q[1], r_l_q[2], k, p, N,
                                                                      s_vector_k, R_vector, m_vector);
                                        H(i_mat, j_mat) += J_ratio * val;
                                    }
                                }
                            }
                        }
                    }
                }
                H_subSubSubspace_list.emplace_back(H);
            }
            H_subSubspace_list.emplace_back(H_subSubSubspace_list);
        }
        H_subspace_list.emplace_back(H_subSubspace_list);
    }
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
    if ((m + 1) < epsilon) {
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


std::vector<int> checkState_parity(const int s, const int k, const int N) {
    int t = s;
    int R = -1;
    int m = -1;
    for (int i = 1; i <= N/2; i++) {
        cycleBits2(t, N); //translate state
        if (t < s) {return {R, m};}
        else if (t == s) {
            if ( k % (int) trunc(N/(2*i)) ) {return {R, m};} // check compatibility with k
            R = i;
            break;
        }
    }
    t = s;

    reflectBits(t, N);
    for (int i = 0; i < R; i++) {
        if (t < s) {
            R = -1;
            return {R, m};
        } else if (t == s) {
            m = i;
            return {R, m};
        }
        cycleBits2(t, N);
    }
    return {R, m};
}

std::vector<int> representative_parity(const int s, const int N) {
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
    return {r, l, q};
}
