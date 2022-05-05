#include "main.h"

#include <chrono>
#include <utility>
#include <thread>
#include <mutex>

using std::vector;
using std::complex;
using std::list;
using std::pow;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;

//#define saveErgs
//#define saveExcitationErgs
//#define saveSpecificHeat
//#define parallelDiag
//#define susceptibility
#define spinTest2

int main(int argc, char* argv[]) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
#ifdef saveErgs
    int N; double j_ratio; std::string path, method;
    if (argc > 4) {
        N = atoi(argv[1]);
        j_ratio = atof(argv[2]);
        path =argv[3];
        method = argv[4];
    } else {
        N = 10;
        j_ratio = 1.0;
        path = "ergs.txt";
        method = "momentum";
    }

    if (method == "naive") {
        MatrixXd H = naiveHamiltonian(j_ratio, N);
        Eigen::VectorXd erg = H.eigenvalues().real();
        std::sort(erg.begin(), erg.end());
        //saveEnergies(erg , path);
    } else if (method == "magnetization") {
        list<MatrixXd> H = magnetizationHamiltonian(j_ratio, N);
        vector<double> erg = getEnergiesFromBlocks(H, N);
        saveEnergies(erg, path);
    } else {
        list<list<MatrixXcd>> H = momentumHamiltonian(j_ratio, N);
        vector<double> erg = getEnergiesFromBlocks(H, N);
        saveEnergies(erg, path);
    }

#endif
#ifdef saveExcitationErgs
    Eigen::VectorXd J_ratios = Eigen::VectorXd::LinSpaced(200, 0, 2.5);
    int N = 10;
    vector<double> diffs;
    for (double J_ratio : J_ratios) {
        list<list<MatrixXcd>> H = momentumHamiltonian(J_ratio, N);
        vector<double> erg = getEnergiesFromBlocks(H, N);
        diffs.emplace_back( abs(erg[0]-erg[1]) );
    }
    list<std::pair<double, double>> out;
    for (int i = 0; i < diffs.size(); i++) {
        out.emplace_back( std::pair<double, double>(J_ratios[i], diffs[i]) );
    }
    std::ofstream ergFile;
    ergFile.open("/home/mmaschke/BA_Code/Data/excitationEnergies.txt");
    for (std::pair<double, double> p : out) {
        ergFile << p.first << " " << p.second << "\n";
    }
    ergFile.close();
#endif
#ifdef saveSpecificHeat
    int dataPointNum = 200;
    Eigen::VectorXd Ts = Eigen::VectorXd::LinSpaced(dataPointNum, 0.001, 2.5);
    double J_ratio = 2;
    int N = 18;
    vector<double> C(dataPointNum);
    int i = 0;
    list<list<MatrixXcd>> H = momentumHamiltonian(J_ratio, N);
    vector<double> erg = getMomentumErgsThreaded(H, N);
    for (double beta : Ts) {
        C[i] = specificHeat(erg, beta, false)/N ;
        i++;
    }
    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(Ts[j], C[j]) );
    }
    std::ofstream ergFile;
    ergFile.open("/home/mmaschke/BA_Code/Data/specHeatJ2N18.txt");
    for (std::pair<double, double> p : out) {
        ergFile << p.first << " " << p.second << "\n";
    }
    ergFile.close();
#endif
#ifdef parallelDiag
    int dataPointNum = 100;
    Eigen::VectorXd J_ratios = Eigen::VectorXd::LinSpaced(dataPointNum, 0, 2.5);
    vector<double> J_ratios_stl(J_ratios.data(), J_ratios.data() + J_ratios.size());
    int N = 14;
    vector<vector<double>> energies_for_different_betas = diagonalizeThreaded(J_ratios_stl, N);
#endif
#ifdef susceptibility
    int N = 8;
    double j_ratio = 0;
    int dataPointNum = 200;
    Eigen::VectorXd Ts = Eigen::VectorXd::LinSpaced(dataPointNum, 0, 2.5);
    MatrixXd S_2 = spinOperator_sq(N);

    MatrixXd H = naiveHamiltonian(j_ratio, N);
    Eigen::ComplexEigenSolver<MatrixXd> sol(H);
    Eigen::VectorXd erg = sol.eigenvalues().real();
    vector<double> ergs_stl(erg.data(), erg.data() + erg.size());
    const Eigen::MatrixXcd & U = sol.eigenvectors();

    vector<double> susceptibilities(dataPointNum);
    for (int i = 0; i < dataPointNum; i++) {
        susceptibilities[i] =  susceptibilityOld(ergs_stl, Ts[i], false, U, S_2) / (double) N ;
    }

    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(Ts[j], susceptibilities[j]) );
    }
    std::ofstream ergFile;
    ergFile.open("/home/mmaschke/BA_Code/Data/suscN8J0beta.txt");
    for (std::pair<double, double> p : out) {
        ergFile << p.first << " " << p.second << "\n";
    }
    ergFile.close();

    //printMatrix(U);
#endif
#ifdef spinTest2
    int N = 8;
    double j_ratio = 0;
    int dataPointNum = 200;
    Eigen::VectorXd Ts = Eigen::VectorXd::LinSpaced(dataPointNum, 0.001, 2.5);
    vector<int> states = getStates_m(N, N/2);
    MatrixXd S_2 = spinOperator_sq(states, N);

    MatrixXd H_m0 = getMagnetizationBlock(j_ratio, 0, N);
    Eigen::ComplexEigenSolver<MatrixXd> sol(H_m0);
    Eigen::VectorXd erg = sol.eigenvalues().real();
    vector<double> ergs_stl(erg.data(), erg.data() + erg.size());
    const Eigen::MatrixXcd & U = sol.eigenvectors();

    vector<double> susceptibilities(dataPointNum);
    for (int i = 0; i < dataPointNum; i++) {
        susceptibilities[i] =  susceptibility(ergs_stl, Ts[i], false, U, S_2) / (double) N ;
    }


    list<std::pair<double, double>> out;
    for (int j = 0; j < dataPointNum; j++) {
        out.emplace_back( std::pair<double, double>(Ts[j], susceptibilities[j]) );
    }
    std::ofstream ergFile;
    ergFile.open("/home/mmaschke/BA_Code/Data/suscNewMethod.txt");
    for (std::pair<double, double> p : out) {
        ergFile << p.first << " " << p.second << "\n";
    }
    ergFile.close();


#endif
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count()
              << "[ms]" << std::endl;
    return 0;
}

vector<double> getMomentumErgsThreaded(const list<list<MatrixXcd>> & H_list, int N) {
    vector<list<MatrixXcd>> H_vector(H_list.begin(), H_list.end());
    vector<double> ergs;
#pragma omp parallel for default(none) shared(ergs, H_vector, N, std::cout) num_threads(16)
    for (int i = 0; i < H_vector.size(); i++) {
        vector<double> blockErgs = getEnergiesFromBlocks(H_vector[i]);
        writeThreadSafe(ergs, blockErgs);
        std::cout << "1 done" << "\n";
    }
    std::sort(ergs.begin(), ergs.end());
    return ergs;
}

vector<vector<double>> diagonalizeThreaded(const vector<double> & J_ratios, int N) {
    const int num = J_ratios.size();
    vector<vector<double>> v(num);
#pragma omp parallel for default(none) shared(v, J_ratios, N, std::cout) num_threads(16)
    for (int i = 0; i < num; i++) {
        list<list<MatrixXcd>> H = momentumHamiltonian(J_ratios[i], N);
        vector<double> erg = getEnergiesFromBlocks(H, N);
        writeThreadSafe(v, erg);
        std::cout << "1 done" << "\n";
    }
    return v;
}

void writeThreadSafe (vector<vector<double>> & writeTo, const vector<double> & writeFrom) {
    static std::mutex mu;
    std::lock_guard<std::mutex> lock(mu);
    writeTo.emplace_back(writeFrom);
}

void writeThreadSafe (vector<double> & writeTo, const vector<double> & writeFrom) {
    static std::mutex mu;
    std::lock_guard<std::mutex> lock(mu);
    for (double d : writeFrom) {
        writeTo.emplace_back(d);
    }
}

