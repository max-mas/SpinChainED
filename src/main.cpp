#include "main.h"

using std::vector;

//#define saveExcitationErgs
//#define saveGroundStateErgsPerSpin
//#define saveSpecificHeat
//#define saveSpecificHeatForJ
//#define saveSusceptibility
//#define saveSusceptibilityForJ
//#define saveDispersion

int main(int argc, char* argv[]) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    if (argc < 10) {
        std::cout <<
            "Correct usage: ./SpinChainEd nMin nMax dataPointNum path_J_ratios path_Ts isBeta startP endP flags.\n"
            "Flag order: ExcErgs, GrdState, SpecHeat, SpecHeatJ, Susc, SuscJ, Disp.\n";
        return 69;
    }

    int nMin = atoi(argv[1]);
    int nMax = atoi(argv[2]);
    int dataPointNum = atoi(argv[3]);
    vector<double> J_ratios = readDoubleVectorFromFile(argv[4]);
    vector<double> Ts = readDoubleVectorFromFile(argv[5]);
    bool isBeta = atoi(argv[6]);
    double start = atof(argv[7]);
    double endP = atof(argv[8]);
    std::string flags_str = argv[9];
    vector<bool> flags(7, false);

    for (int i = 0; i < flags_str.length(); i++) {
        char flag = flags_str[i] - '0';
        flags[i] = (int) flag;
    }

    // Excitation energies
    if (flags[0]) {
        for (int N = nMin; N <= nMax; N += 2) {
            std::string path = "/home/mmaschke/BA_Code/Data/ExcitationErgs/ExcErgs" + std::to_string(N) + ".txt";
            saveExcitationErgsForVaryingJ(N, dataPointNum, start, endP, path);
            std::cout << std::string("N") + std::to_string(N) << std::endl;
        }
    }

    // Ground state energies
    if (flags[1]) {
        for (int N = nMin; N <= nMax; N+= 2) {
            std::string path = "/home/mmaschke/BA_Code/Data/GroundStateErgs/GSErgs" + std::to_string(N) + ".txt";
            saveGroundStateErgPerSpinForVaryingJ(N, dataPointNum, start, endP, path);
            std::cout << std::string("N") + std::to_string(N) << std::endl;
        }
    }

    // Specific heats (T)
    if (flags[2]) {
        std::cout << "Specific heats vor varying temps:\n";
        for (double J_ratio: J_ratios) {
            for (int N = nMin; N <= nMax; N += 2) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path =
                        "/home/mmaschke/BA_Code/Data/SpecificHeats/SpecHeatN" + std::to_string(N) + std::string("J") +
                        j + ".txt";
                saveSpecificHeatsForVaryingTemp(N, dataPointNum, J_ratio, start, endP, isBeta, path);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j << std::endl;
            }
        }
    }

    // Specific heats (J)
    if (flags[3]) {
        std::cout << "Specific heats vor varying J:\n";
        for (double T: Ts) {
            for (int N = nMin; N <= nMax; N += 2) {
                std::string b = std::to_string(T);
                std::replace(b.begin(), b.end(), '.', '_');
                std::string path = "/home/mmaschke/BA_Code/Data/SpecificHeatsForJ/SpecHeatN" + std::to_string(N) +
                                   std::string("T") +
                                   b + ".txt";
                saveSpecificHeatsForVaryingJ(N, dataPointNum, T, start, endP, isBeta, path);
                std::cout << std::string("N") + std::to_string(N) + std::string("T") + b << std::endl;
            }
        }
    }

    // Susceptibilities (T)
    if (flags[4]) {
        std::cout << "Susceptibilities vor varying temps:\n";
        for (double J_ratio: J_ratios) {
            for (int N = nMin; N <= nMax; N += 2) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path =
                        "/home/mmaschke/BA_Code/Data/Susceptibilities/SuscN" + std::to_string(N) + std::string("J") +
                        j + ".txt";
                saveSusceptibilitesForVaryingTemp(N, dataPointNum, J_ratio, start, endP, isBeta, path);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j << std::endl;
            }
        }
    }

    // Susceptibilities (J)
    if (flags[5]) {
        std::cout << "Susceptibilities vor varying J:\n";
        for (double T : Ts) {
            for (int N = nMin; N <= nMax; N+= 2) {
                std::string b = std::to_string(T);
                std::replace(b.begin(), b.end(), '.', '_');
                std::string path = "/home/mmaschke/BA_Code/Data/SusceptibilitiesForJ/SuscN" + std::to_string(N)+ std::string("T") +
                                   b + ".txt";
                saveSusceptibilitiesForVaryingJ(N, dataPointNum, T, start, endP, isBeta, path);
                std::cout << std::string("N") + std::to_string(N)+ std::string("T") + b << std::endl;
            }
        }
    }

    // Dispersion
    if (flags[6]) {
        for (double J_ratio : J_ratios) {
            for (int N = nMin; N <= nMax; N+= 2) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');

                std::string path = "/home/mmaschke/BA_Code/Data/Dispersion/DispN" + std::to_string(N)
                                   + std::string("J") + j + ".txt";
                saveEnergyDispersionWithMag(N, J_ratio, path);
                std::cout << std::string("N") + std::to_string(N)+ std::string("J") + j << std::endl;
            }
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count()
              << "[ms]" << std::endl;
    return 0;
}
