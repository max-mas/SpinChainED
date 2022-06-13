#include "main.h"

using std::vector;

//#define saveExcitationErgs
//#define saveGroundStateErgsPerSpin
//#define saveSpecificHeat
//#define saveSpecificHeatForJ
//#define saveSusceptibility
//#define saveSusceptibilityForJ
//#define saveDispersion
#define QTtestingArea
#define statisticsTest

int main(int argc, char* argv[]) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

#ifndef QTtestingArea
    if (argc < 11) {
        std::cout <<
            "Correct usage: ./SpinChainEd nMin nMax dataPointNum path_J_ratios path_Ts isBeta startP endP flags saveTo_path.\n"
            "Flag order: ExcErgs, GrdState, SpecHeat, SpecHeatJ, Susc, SuscJ, Disp.\n";
        return 69;
    }

    const int cpu_cnt = (int) std::thread::hardware_concurrency() ;// / 2;
    omp_set_num_threads(cpu_cnt);


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

    std::string saveTo_path = argv[10];

    for (int i = 0; i < flags_str.length(); i++) {
        char flag = flags_str[i] - '0';
        flags[i] = (int) flag;
    }

    // Excitation energies
    if (flags[0]) {
        std::cout << "Excitation energies:\n";
        for (int N = nMin; N <= nMax; N += 2) {
            std::string path = saveTo_path + "/out/ExcitationErgs/ExcErgs" + std::to_string(N) + ".txt";
            saveSpinGapForVaryingJ(N, dataPointNum, start, endP, path);
            std::cout << std::string("N") + std::to_string(N) << std::endl;
        }
    }

    // Ground state energies
    if (flags[1]) {
        std::cout << "Ground state energies:\n";
        for (int N = nMin; N <= nMax; N+= 2) {
            std::string path = saveTo_path + "/out/GroundStateErgs/GSErgs" + std::to_string(N) + ".txt";
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
                std::string path = saveTo_path + "/out/SpecificHeats/SpecHeatN" + std::to_string(N) + std::string("J") +
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
                std::string path = saveTo_path + "/out/SpecificHeatsForJ/SpecHeatN" + std::to_string(N) +
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
                std::string path = saveTo_path + "/out/Susceptibilities/SuscN" + std::to_string(N) + std::string("J") +
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
                std::string path = saveTo_path + "/out/SusceptibilitiesForJ/SuscN" + std::to_string(N)+ std::string("T") +
                                   b + ".txt";
                saveSusceptibilitiesForVaryingJ(N, dataPointNum, T, start, endP, isBeta, path);
                std::cout << std::string("N") + std::to_string(N)+ std::string("T") + b << std::endl;
            }
        }
    }

    // Dispersion
    if (flags[6]) {
        for (double J_ratio : J_ratios) {
            std::cout << "Dispersion:\n";
            for (int N = nMin; N <= nMax; N+= 2) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');

                std::string path = saveTo_path + "/out/Dispersion/DispN" + std::to_string(N)
                                   + std::string("J") + j + ".txt";
                saveEnergyDispersionWithMag(N, J_ratio, path);
                std::cout << std::string("N") + std::to_string(N)+ std::string("J") + j << std::endl;
            }
        }
    }
#endif
#ifdef QTtestingArea

    vector<double> J_ratios = {0.1, 0.5, 1, 2};
    vector<int> runNums = {1, 3, 10};
    int nMin = 14;
    int nMax = 12;
    std::string saveTo_path = "D:/Code/C++/spinChainData";
    int dataPointNum = 1000;

    for (double J_ratio: J_ratios) {
        for (int numOfRuns : runNums) {
            for (int N = nMin; N <= nMax; N += 2) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path = saveTo_path + "/out/SpecificHeats_DQT/SpecHeatDQTN" + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                saveSpecificHeatsForVaryingTemp_DQT_avg(N, dataPointNum, J_ratio, 10, path, numOfRuns);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j + std::string("It")
                        + std::to_string(numOfRuns) << std::endl;
            }
        }
    }

#endif
#ifdef statisticsTest
    for (double J_ratio: J_ratios) {
        for (int numOfRuns : runNums) {
            for (int N = 6; N <= 16; N += 2) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string EDPath = saveTo_path + "/out/SpecificHeats/SpecHeatN" + std::to_string(N)
                                   + std::string("J") + j + ".txt";
                std::string DQTPath = saveTo_path + "/out/SpecificHeats_DQT/SpecHeatDQTN" + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                std::string outPath = saveTo_path + std::string("/out/QTErrorStats/SpecHeatDiffs/DiffN") + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                calcAbsDQTError(EDPath, DQTPath, outPath);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j + std::string("It")
                             + std::to_string(numOfRuns) << std::endl;
            }
        }
    }
#endif

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count()
              << "[ms]" << std::endl;
    return 0;
}
