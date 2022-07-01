#include "main.h"

#include <sys/resource.h>
#include <malloc.h>

using std::vector;

//#define saveExcitationErgs
//#define saveGroundStateErgsPerSpin
//#define saveSpecificHeat
//#define saveSpecificHeatForJ
//#define saveSusceptibility
//#define saveSusceptibilityForJ
//#define saveDispersion
#define QTtestingArea
//#define QTDataForFit
//#define statisticsTest
#define EDbenchmark

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
            std::string path = saveTo_path + "/out/ActualExcitationErgs/ExcErgs" + std::to_string(N) + ".txt";
            saveExcitationErgsForVaryingJ(N, dataPointNum, start, endP, path);
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
        for (int N = nMin; N <= nMax; N += 2) {
            for (double J_ratio: J_ratios) {
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
    vector<int> runNums = {2, 5};
    int nMin = 16;
    int nMax = 12;
    std::string saveTo_path = "/home/mmaschke/BA_Code/Data";
    int dataPointNum = 5000;

    /*
    for (int N = nMin; N <= nMax; N += 2) {
#pragma omp parallel for default(none) shared(runNums, saveTo_path, N, dataPointNum, std::cout, J_ratios)
        for (int l = 0; l < J_ratios.size(); l++) {
            double J_ratio = J_ratios[l];
            for (int k = 0; k < runNums.size(); k++) {
                int numOfRuns = runNums[k];
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path = saveTo_path + "/out/SpecificHeats_DQT/SpecHeatDQTN" + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                saveSpecificHeatsForVaryingTemp_DQT_avg(N, dataPointNum, J_ratio, 50, path, numOfRuns);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j + std::string("It")
                             + std::to_string(numOfRuns) << std::endl;
            }
        }
    }*/

    for (int N = nMin; N <= nMax; N += 2) {
        const Eigen::SparseMatrix<std::complex<double>> S2 = spinOp2_momentum_sparse(N);
#pragma omp parallel for default(none) shared(runNums, saveTo_path, N, dataPointNum, std::cout, J_ratios)
        for (int l = 0; l < J_ratios.size(); l++) {
            double J_ratio = J_ratios[l];
            for (int k = 0; k < runNums.size(); k++) {
                int numOfRuns = runNums[k];
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path = saveTo_path + "/out/Susceptibilities_DQT/SuscDQTN" + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                saveSusceptibilityForVaryingTemp_DQT_avg(N, dataPointNum, J_ratio, 50, S2, path, numOfRuns);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j + std::string("It")
                                + std::to_string(numOfRuns) << std::endl;
            }
        }
    }

#endif
#ifdef QTDataForFit
    Eigen::VectorXd Js = Eigen::VectorXd::LinSpaced(50, 0, 2);
    nMin = 6;
    nMax = 12;
    int numOfRuns = 1;

    /*
    for (int N = nMin; N <= nMax; N += 2) {
        for (int i = 1; i <= numOfRuns; i++) {
#pragma omp parallel for default(none) shared(Js, saveTo_path, dataPointNum, std::cout, i, N)
            for (int k = 0; k < Js.size(); k++) {
                double J_ratio = Js[k];
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path = saveTo_path + "/out/SpecificHeats_DQT/forFit/" + std::to_string(i) + "/SpecHeatDQTN" + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(1) + ".txt";
                saveSpecificHeatsForVaryingTemp_DQT_parallel(N, dataPointNum, J_ratio, 50, path);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j + std::string("It")
                             + std::to_string(1) << std::endl;
            }
        }
    }

    for (int N = nMin; N <= nMax; N += 2) {
        const Eigen::SparseMatrix<std::complex<double>> S2 = spinOp2_momentum_sparse(N);
        for (int i = 1; i <= numOfRuns; i++) {
            for (double J_ratio: Js) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path = saveTo_path + "/out/Susceptibilities_DQT/forFit/" + std::to_string(i) + "/SuscDQTN" + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(1) + ".txt";
                saveSusceptibilityForVaryingTemp_DQT_avg(N, dataPointNum, J_ratio, 50, S2, path, 1);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j + std::string("It")
                             + std::to_string(1) << std::endl;
            }
        }
    }

    for (int i = 1; i <= numOfRuns; i++) {
        for (int N = nMin; N <= nMax; N += 2) {
            for (double J_ratio: Js) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path = saveTo_path + "/out/SpecificHeats/forFit/" + std::to_string(i) + "/SpecHeatN"
                        + std::to_string(N) + std::string("J") + j + ".txt";
                saveSpecificHeatsForVaryingTemp(N, dataPointNum, J_ratio, 0, 50, true, path);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j << std::endl;
            }
        }
    }

    for (int i = 1; i <= numOfRuns; i++) {
        for (int N = nMin; N <= nMax; N += 2) {
            for (double J_ratio: Js) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string path = saveTo_path + "/out/Susceptibilities/forFit/" + std::to_string(i) + "/SuscN"
                        + std::to_string(N) + std::string("J") + j + ".txt";
                saveSusceptibilitesForVaryingTemp(N, dataPointNum, J_ratio, 0, 50, true, path);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j << std::endl;
            }
        }
    }*/


#endif
#ifdef statisticsTest
    for (int N = 6; N <= 16; N += 2) {
        for (double J_ratio: J_ratios) {
            for (int numOfRuns : runNums) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string EDPath = saveTo_path + "/out/SpecificHeats/SpecHeatN" + std::to_string(N)
                                     + std::string("J") + j + ".txt";
                std::string DQTPath = saveTo_path + "/out/SpecificHeats_DQT/SpecHeatDQTN" + std::to_string(N)
                                      + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                std::string outPath = saveTo_path + std::string("/out/QTErrorStats/SpecHeatDiffs/DiffN") + std::to_string(N)
                                      + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                calcRelDQTError(EDPath, DQTPath, outPath);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j + std::string("It")
                             + std::to_string(numOfRuns) << std::endl;
            }
        }
    }

    for (int N = 6; N <= 16; N += 2) {
        for (double J_ratio: J_ratios) {
            for (int numOfRuns : runNums) {
                std::string j = std::to_string(J_ratio);
                std::replace(j.begin(), j.end(), '.', '_');
                std::string EDPath = saveTo_path + "/out/Susceptibilities/SuscN" + std::to_string(N)
                                   + std::string("J") + j + ".txt";
                std::string DQTPath = saveTo_path + "/out/Susceptibilities_DQT/SuscDQTN" + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                std::string outPath = saveTo_path + std::string("/out/QTErrorStats/SuscDiffs/DiffN") + std::to_string(N)
                                   + std::string("J") + j + std::string("It") + std::to_string(numOfRuns) + ".txt";
                calcRelDQTError(EDPath, DQTPath, outPath);
                std::cout << std::string("N") + std::to_string(N) + std::string("J") + j + std::string("It")
                             + std::to_string(numOfRuns) << std::endl;
            }
        }
    }
#endif
#ifdef EDbenchmark
    int minN = 6;
    int maxN = 16;
    std::chrono::steady_clock::time_point start;
    std::chrono::steady_clock::time_point finish;
    /*
    std::cout << "Momentum:" << std::endl;
    for (int N = minN; N <= maxN; N+=2) {
        start = std::chrono::steady_clock::now();
        std::list<std::list<Eigen::MatrixXcd>> H = momentumHamiltonian(0, N, 0, N);
        vector<double> eig = getEnergiesFromBlocks(H, true);
        finish = std::chrono::steady_clock::now();
        std::ifstream statFile("/proc/self/stat");
        std::string statLine;
        std::getline(statFile, statLine);
        std::istringstream iss(statLine);
        std::string entry;
        long long memUsage;
        for (int i = 1; i <= 24; i++) {
            std::getline(iss, entry, ' ');
            if (i == 24) {
                memUsage = stoi(entry);
            }
        }

        long time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count();
        std::cout << "N = " << N << ", time = " << time_ns/1e9 << " s, ram = " << 4096*memUsage/1e9 << " GB" << std::endl;
    }*/

    std::cout << "DQT Specific Heat, 5000 Data Points:" << std::endl;
    for (int N = minN; N <= maxN; N+=2) {
        saveSpecificHeatsForVaryingTemp_DQT_avg(N, 5000, 0, 50, "", 1, true);
    }

#endif

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << "[ms]" << std::endl;
    return 0;
}
