#include "main.h"
#include <chrono>

using std::vector;
using std::complex;
using std::list;
using std::pow;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;

int main(int argc, char* argv[]) {

    int N = 6;
    double j_ratio = 2;

    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


    MatrixXd H = naiveHamiltonian(j_ratio, N);
    //printMatrix(H);
    Eigen::VectorXd erg = H.eigenvalues().real();
    std::sort(erg.begin(), erg.end());
    printEnergies(erg);
    //printMatrix(H);

    std::cout << std::endl;

    list<MatrixXd> H2 = magnetizationHamiltonian(j_ratio, N);
    vector<double> erg2 = getEnergiesFromBlocks(H2, N);
    printEnergies(erg2);

    std::cout << std::endl;

    list<list<MatrixXcd>> H3 = momentumHamiltonian(j_ratio, N);
    vector<double> erg3 = getEnergiesFromBlocks(H3, N);
    printEnergies(erg3);
    //list<MatrixXcd> H3_1 = blkdiag(H3);
    //MatrixXcd H3_2 = blkdiag(H3_1, pow(2, N));
    //printMatrix(H3_2);
    //MatrixXcd H3_2_tr = H3_2.transpose();
    //printMatrix(H3_2_tr);

    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;

    return 0;
}

