#include "main.h"

using std::vector;

//#define saveErgs
//#define saveExcitationErgs
//#define saveSpecificHeat
//#define saveSpecificHeatForJ
//#define saveSusceptibility
//#define saveSusceptibilityForJ
#define saveDispersion

int main(int argc, char* argv[]) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    int dataPointNum = 200;
    int nMax = 12;
    vector<double> J_ratios = {0.1, 1, 2, 5};
    vector<double> Ts = {0.2, 0.5, 1};
    bool isBeta = false;
    double start = 0;
    double endP = 3;
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
    for (int N = 6; N <= nMax; N+= 2) {
        std::string path = "/home/mmaschke/BA_Code/Data/ExcitationErgs/ExcErgs" + std::to_string(N) + ".txt";
        saveExcitationErgsForVaryingJ(N, dataPointNum, start, endP, path);
        std::cout << std::string("N") + std::to_string(N) << std::endl;
    }
#endif
#ifdef saveSpecificHeat
    for (double J_ratio : J_ratios) {
        for (int N = 6; N <= nMax; N+= 2) {
            std::string j = std::to_string(J_ratio);
            std::replace(j.begin(), j.end(), '.', '_');
            std::string path = "/home/mmaschke/BA_Code/Data/SpecificHeats/SpecHeatN" + std::to_string(N)+ std::string("J") +
                    j + ".txt";
            saveSpecificHeatsForVaryingTemp(N, dataPointNum, J_ratio, start, endP, isBeta, path);
            std::cout << std::string("N") + std::to_string(N)+ std::string("J") + j << std::endl;
        }
    }
#endif
#ifdef saveSpecificHeatForJ
    for (double T : Ts) {
        for (int N = 6; N <= nMax; N+= 2) {
            std::string b = std::to_string(T);
            std::replace(b.begin(), b.end(), '.', '_');
            std::string path = "/home/mmaschke/BA_Code/Data/SpecificHeatsForJ/SpecHeatN" + std::to_string(N)+ std::string("T") +
                               b + ".txt";
            saveSpecificHeatForVaryingJ(N, dataPointNum, T, start, endP, isBeta, path);
            std::cout << std::string("N") + std::to_string(N)+ std::string("T") + b << std::endl;
        }
    }
#endif
#ifdef saveSusceptibility
    for (double J_ratio : J_ratios) {
        for (int N = 6; N <= nMax; N+= 2) {
            std::string j = std::to_string(J_ratio);
            std::replace(j.begin(), j.end(), '.', '_');
            std::string path = "/home/mmaschke/BA_Code/Data/Susceptibilities/SuscN" + std::to_string(N)+ std::string("J") +
                               j + ".txt";
            saveSusceptibilitesForVaryingTemp(N, dataPointNum, J_ratio, start, endP, isBeta, path);
            std::cout << std::string("N") + std::to_string(N)+ std::string("J") + j << std::endl;
        }
    }
#endif
#ifdef saveSusceptibilityForJ
    for (double T : Ts) {
        for (int N = 6; N <= 12; N+= 2) {
            std::string b = std::to_string(T);
            std::replace(b.begin(), b.end(), '.', '_');
            std::string path = "/home/mmaschke/BA_Code/Data/SusceptibilitiesForJ/SuscN" + std::to_string(N)+ std::string("T") +
                               b + ".txt";
            saveSusceptibilitiesForVaryingJ(N, dataPointNum, T, start, endP, isBeta, path);
            std::cout << std::string("N") + std::to_string(N)+ std::string("T") + b << std::endl;
        }
    }
#endif
#ifdef saveDispersion
    for (double J_ratio : J_ratios) {
        for (int N = 6; N <= nMax; N+= 2) {
            std::string j = std::to_string(J_ratio);
            std::replace(j.begin(), j.end(), '.', '_');

            std::string path = "/home/mmaschke/BA_Code/Data/Dispersion/DispN" + std::to_string(N)
                    + std::string("J") + j + ".txt";
            saveEnergyDispersion(N, J_ratio, path);
            std::cout << std::string("N") + std::to_string(N)+ std::string("J") + j << std::endl;
        }
    }
#endif
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count()
              << "[ms]" << std::endl;
    return 0;
}
