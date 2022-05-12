//
// This cpp file contains functions used to manipulate integers in the way needed for the chosen numerical
// representation of spin states, as well as functions used for console output.
//
#include <iostream>

#include "helperFunctions.h"

using std::vector;
using std::complex;

using Eigen::MatrixXcd;
using Eigen::Dynamic;


// Gets bit in integer s at index i.
bool getBit(const int s, const int i) {
    if (s & (1 << i)) {
        return true;
    } else {
        return false;
    }
}

// Flips bit in int s at index i using bitwise XOR.
void flipBit(int &s, const int i) {
    s ^= (1 << i);
}

// Sets bit in int s at index i using bitwise AND.
void setBit(int & s, const int i, const bool val) {
    if (val) {
        s |= (1 << i);
    } else {
        s &= ~(1 << i);
    }
}

// Counts set (=1) bits in int s using bitwise AND.
int bitSum(const int s) {
    int count = 0;
    for (int i = 0; i <= 31; i++) {
        if (getBit(s, i)) {
            count++;
        }
    }
    return count;
}

// Returns index of first occurrence of int s in s_list. Returns -1 if not found.
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

// Calls cycleBits twice.
void cycleBits2(int &s, const int N) {
    cycleBits(s, N);
    cycleBits(s, N);
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

// Reflects bits about the center of the chain.
void reflectBits(int & s, int N) {
    int t = 0;
    for (int i = 0; i < N; i++) {
        setBit(t, i, getBit(s, N - 1 - i));
    }
    s = t;
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

std::vector<int> representative_parity(const int s, const int N) {
    int r = s;
    int t = s;
    int l = 0;

    for (int i = 1; i < N/2; i++) {
        cycleBits2(t, N);
        if (t < r) {r = t; l = i;}
    }

    t = s;
    reflectBits(t, N);
    int q = 0;
    for (int i = 1; i < N/2; i++) {
        cycleBits2(t, N);
        if (t < r) {
            r = t;
            l = i;
            q = 1;
        }
    }
    return {r, l, q};
}

// Prints real-valued matrix to console.
void printMatrix(const Eigen::MatrixXd & M) {
    for (int i = 0; i < M.rows(); i++ ) {
        for (int j = 0; j < M.cols(); j++) {
            if (M(i,j) < epsilon) {
                std::cout << 0.0 << "\t\t\t" ;
            } else {
                std::cout << M(i, j) << "\t\t\t" ;
            }

        }
        std::cout << std::endl;
    }
}

// Prints complex-valued matrix to console.
void printMatrix(const Eigen::MatrixXcd & M) {
    for (int i = 0; i < M.rows(); i++ ) {
        for (int j = 0; j < M.cols(); j++) {
            if (abs(M(i,j)) < epsilon) {
                std::cout << 0.0 << "\t\t\t" ;
            } else {
                std::cout << M(i, j) << "\t\t\t" ;
            }
        }
        std::cout << std::endl;
    }
}

// Prints doubles to console from Eigen::Vector.
void printEnergies(const Eigen::VectorXd & v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << v(i) << std::endl;
    }
}

// Prints doubles to console from std::vector.
void printEnergies(const vector<double> & v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << v[i] << std::endl;
    }
}