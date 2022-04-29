//
// Created by mmaschke on 26/04/22.
//
#include <iostream>

#include "helperFunctions.h"

using std::vector;
using std::complex;

using Eigen::MatrixXcd;
using Eigen::Dynamic;

// Helper functions

// Get bit in integer s at index i.
bool getBit(const int s, const int i) {
    if (s & (1 << i)) {
        return true;
    } else {
        return false;
    }
}

// Flip bit in int s at index i using bitwise XOR.
void flipBit(int &s, const int i) {
    s ^= (1 << i);
}

// Set bit in int s at index i using bitwise AND.
void setBit(int & s, const int i, const bool val) {
    if (val) {
        s |= (1 << i);
    } else {
        s &= ~(1 << i);
    }
}

// Count set (=1) bits in int s using bitwise AND.
int bitSum(const int s) {
    int count = 0;
    for (int i = 0; i <= 31; i++) {
        if (getBit(s, i)) {
            count++;
        }
    }
    return count;
}

// Return index of first occurrence of int s in s_list. Return -1 if not found.
int findState(const std::vector<int> &s_list, const int s) {
    for (int i = 0; i < s_list.size(); i++) {
        if ( s_list[i] == s ) {return i;}
    }
    return -1;
}

// Cyclic left bit shift by 1 for integers in which only the first N bits must be set! Otherwise, undefined behavior.
void cycleBits(int &s, const int N) {
    int i = N - 1;
    // Find first set bit from the left.
    for (int discard = 0; i >= 0; i--) {
        if ( s & (1 << i) ) {
            break;
        }
    }
    // Special case if s is 1 or 0;
    if ( (i == 0) && s == 1 ) { s = 2; return; }
    if ( (i == 0) && s == 0 ) { s = 0; return; }

    // Shift bits cyclically.
    s = s << 1;
    if (i != N - 1) {
        return;
    } else {
        setBit(s, 0, true);
        setBit(s, i+1, false);
    }
}

// Checks if state int s is smallest among its translations and compatible with momentum k.
int checkState(const int s, const int k, const int N) {
    int t = s;
    for (int i = 1; i <= N; i++) {
        cycleBits(t, N); //translate state
        if (t < s) {return -1;}
        else if (t == s) {
            if ( k % (int) std::trunc((double )N/i) ) {return -1;} // check compatibility with k
            return i;
        }
    }
    return -1;
}

// jfc
vector<int> checkState_parity(int s, int k, int N) {
    vector<int> ret = {-1, -1};
    // ret[0] = R, ret[1] = m
    int t = s;
    for (int i = 1; i<= N; i++) {
        cycleBits(t, N);
        if (t < s) { return ret;}
        else if (t == s) {
            if(k % (int) std::trunc((double )N/i)) {return ret;};
            ret[0] = i;
        }
    }
    t = s;
    reflectBits(t, N);
    for (int i = 0; i < ret[0]; i++) {
        if (t < s) {ret[0] = -1; return ret;
        } else if (t == s) {
            ret[1] = i;
        }
        cycleBits(t, N);
    }
    return ret;
}

void reflectBits(int & s, int N) {
    int t;
    for (int i = 0; i < N; i++) {
        setBit(t, i, getBit(s, N - 1 - i));
    }
    s = t;
}

// Generate representative for given state int s. 1st return: representative. 2nd return: Needed number of translations.
std::vector<int> representative(const int s, const int N) {
    int r = s;
    int t = s;
    int l = 0;

    for (int i = 1; i < N; i++) {
        cycleBits(t, N);
        if (t < r) {r = t; l = i;}
    }
    return {r, l};
}

std::vector<int> representative_parity(int s, int N) {
    int r = s;
    int t = s;
    int l = 0;
    int q = 0;

    reflectBits(t, N);
    for (int i = 1; i < N; i++) {
        cycleBits(t, N);
        if (t < r) {r = t; l = i; q = 1;}
    }
    return {r, l, q};
}

double E_z(int s, int N) {
    double E = 0.0;
    for (int i = 0; i < N; i++) {
        int j = (i + 1) % N;
        if (getBit(s, i) == getBit(s, j)) {
            E += 0.25;
        } else {
            E += -0.25;
        }
    }
    return E;
}

double N_sigma(int s, int N, int R, int k, int sigma, int m, int p) {
    double k_actual = 2.0 * M_PI * (double) k / (double) N;
    int g_k;
    if (k_actual > epsilon && M_PI - k_actual > epsilon) {
        g_k = 1;
    } else {
        g_k = 2;
    }
    int s_transl_refl = s;
    reflectBits(s_transl_refl, N);
    double num;
    for (int i = 0; i < N; i++) {
        cycleBits(s_transl_refl, N);
        if (s_transl_refl == s) {
            num = 1.0;
            break;
        }
    }
    if (num != 1.0) {num = 1.0 + (double) sigma * (double) p * cos(k_actual * (double) m);};
    double N_ret = (double) pow(N, 2) * (double) g_k / (double) R * num;
    return N_ret;
}

double hElement(int a, int b, int l, int q, int k, int p, int sigma) {
    return 0;
}

void printMatrix(const Eigen::MatrixXd & M) {
    for (int i = 0; i < M.rows(); i++ ) {
        for (int j = 0; j < M.cols(); j++) {
            std::cout << M(i, j) << "\t" ;
        }
        std::cout << std::endl;
    }
}

void printEnergies(const Eigen::VectorXd & v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << v(i) << std::endl;
    }
}

void printEnergies(const vector<double> & v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << v[i] << std::endl;
    }
}