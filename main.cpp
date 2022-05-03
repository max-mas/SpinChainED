#include "main.h"
#include <chrono>
#include <utility>

using std::vector;
using std::complex;
using std::list;
using std::pow;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::Dynamic;

//#define saveErgs
#define saveExcitationErgs
//#define saveSpecificHeat

int main(int argc, char* argv[]) {

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
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

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

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count()
            << "[ms]" << std::endl;
#endif
#ifdef saveExcitationErgs
    Eigen::VectorXd J_ratios = Eigen::VectorXd::LinSpaced(1000, 0, 5);
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
    ergFile.open("excitationEnergies.txt");
    for (std::pair<double, double> p : out) {
        ergFile << p.first << "\t" << p.second << "\n";
    }
    ergFile.close();

#endif
#ifdef saveSpecificHeat
    list<list<MatrixXcd>> H = momentumHamiltonian(2, 12);
    vector<double> erg = getEnergiesFromBlocks(H, 12);
    std::cout << specificHeat(erg, 1)/12 << "\n";
    return 0;
#endif
}

