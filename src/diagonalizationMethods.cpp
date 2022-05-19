/**
 * This cpp file contains methods used to get Eigenvalues from lists of self-adjoint H-blocks.
*/

#include "diagonalizationMethods.h"

using std::vector;
using std::list;
using std::complex;

using Eigen::MatrixXcd;
using Eigen::MatrixXd;

// Diagonalization methods for a number of cases.
vector<double> getEnergiesFromBlocks(const list<list<list<MatrixXd>>> & H_list, bool sort) {
    vector<double> energies;
    for (const list<list<MatrixXd>> & H_sublist : H_list) {
        for (const list<MatrixXd> & H_subSubList : H_sublist) {
            for (const MatrixXd & mat : H_subSubList) {
                Eigen::SelfAdjointEigenSolver<MatrixXd> sol;
                if (mat.cols() == 0) {
                    continue;
                }
                sol.compute(mat);
                Eigen::VectorXcd blockEnergies = sol.eigenvalues();
                std::for_each(blockEnergies.begin(), blockEnergies.end(),
                        [&](complex<double> & d){energies.emplace_back(d.real());});
            }
        }
    }
    if (sort) std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<MatrixXcd> & H_list, bool sort) {
    vector<double> energies;
    for (const MatrixXcd & mat : H_list) {
        Eigen::SelfAdjointEigenSolver<MatrixXcd> sol;
        if (mat.cols() == 0) {
            continue;
        }
        sol.compute(mat);
        Eigen::VectorXcd blockEnergies = sol.eigenvalues();
        std::for_each(blockEnergies.begin(), blockEnergies.end(),
                      [&](complex<double> & d){energies.emplace_back( d.real() );});
    }
    if (sort) std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<MatrixXd> & H_list, bool sort) {
    vector<double> energies;
    for (const MatrixXd & mat : H_list) {
        Eigen::SelfAdjointEigenSolver<MatrixXd> sol;
        if (mat.cols() == 0) {
            continue;
        }
        sol.compute(mat);
        Eigen::VectorXcd blockEnergies = sol.eigenvalues();
        std::for_each(blockEnergies.begin(), blockEnergies.end(),
                      [&](complex<double> & d){energies.emplace_back( d.real() );});
    }
    if (sort) std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<list<MatrixXcd>> & H_list, bool sort) {
    vector<double> energies;
    for (const list<MatrixXcd> & subList : H_list) {
        for (const MatrixXcd & mat : subList) {
            Eigen::SelfAdjointEigenSolver<MatrixXcd> sol;
            if (mat.cols() == 0) {
                continue;
            }
            sol.compute(mat);
            Eigen::VectorXcd blockEnergies = sol.eigenvalues();
            std::for_each(blockEnergies.begin(), blockEnergies.end(),
                          [&](complex<double> &d) {
                              energies.emplace_back(d.real()); });
        }
    }
    if (sort) std::sort(energies.begin(), energies.end());
    return energies;
}

vector<double> getEnergiesFromBlocks(const list<list<MatrixXd>> & H_list, bool sort) {
    vector<double> energies;
    for (const list<MatrixXd> & subList : H_list) {
        for (const MatrixXd & mat : subList) {
            Eigen::SelfAdjointEigenSolver<MatrixXd> sol;
            if (mat.cols() == 0) {
                continue;
            }
            sol.compute(mat);
            Eigen::VectorXcd blockEnergies = sol.eigenvalues();
            std::for_each(blockEnergies.begin(), blockEnergies.end(),
                          [&](complex<double> &d) { energies.emplace_back(d.real()); });
        }
    }
    if (sort) std::sort(energies.begin(), energies.end());
    return energies;
}

// Returns energies sorted by m and k, used for generation of dispersion plots.
vector<vector<vector<double>>> getEnergiesFromBlocksByK(const list<list<MatrixXcd>> & H_list) {
    vector<vector<vector<double>>> energies;

    for (const list<MatrixXcd> & subList : H_list) {
        vector<vector<double>> m_ergs;
        for (const MatrixXcd & mat : subList) {
            Eigen::SelfAdjointEigenSolver<MatrixXcd> sol;
            if (mat.cols() == 0) {
                continue;
            }
            sol.compute(mat);
            Eigen::VectorXd blockEnergies = sol.eigenvalues().real();
            vector<double> k_ergs(blockEnergies.begin(), blockEnergies.end());
            m_ergs.emplace_back(k_ergs);
        }
        energies.emplace_back(m_ergs);
    }
    return energies;
}

// Threaded version of getEnergiesFromBlocks for the momentum state ansatz.
vector<double> getMomentumErgsThreaded(const list<list<MatrixXcd>> & H_list, int N, bool sort) {
    vector<list<MatrixXcd>> H_vector(H_list.begin(), H_list.end());
    vector<double> ergs;
#pragma omp parallel for default(none) shared(ergs, H_vector, N, std::cout, sort)
    for (int i = 0; i < H_vector.size(); i++) {
        vector<double> blockErgs = getEnergiesFromBlocks(H_vector[i], sort);
        writeThreadSafe(ergs, blockErgs);
        //std::cout << "1 done" << "\n";
    }
    if (sort) std::sort(ergs.begin(), ergs.end());
    return ergs;
}

// Threaded version of getEnergiesFromBlocks for the semi-momentum/parity state ansatz.
vector<double> getParityErgsThreaded(const list<list<list<MatrixXd>>> & H_list, int N, bool sort) {
    vector<list<list<MatrixXd>>> H_vector(H_list.begin(), H_list.end());
    vector<double> ergs;
#pragma omp parallel for default(none) shared(ergs, H_vector, N, std::cout, sort)
    for (int i = 0; i < H_vector.size(); i++) {
        vector<double> blockErgs = getEnergiesFromBlocks(H_vector[i], sort);
        writeThreadSafe(ergs, blockErgs);
        //std::cout << "1 done" << "\n";
    }
    if (sort) std::sort(ergs.begin(), ergs.end());
    return ergs;
}

// Returns all eigenvalues for a given vector of values of J1/J2. Threaded.
vector<vector<double>> diagonalizeThreaded(const vector<double> & J_ratios, int N, bool sort) {
    vector<vector<double>> v(J_ratios.size());
#pragma omp parallel for default(none) shared(v, J_ratios, N, std::cout, sort)
    for (int i = 0; i < J_ratios.size(); i++) {
        list<list<MatrixXcd>> H = momentumHamiltonian(J_ratios[i], N, 0, N);
        vector<double> erg = getEnergiesFromBlocks(H, sort);
        writeThreadSafe(v, erg);
        std::cout << "1 done" << "\n";
    }
    return v;
}

// Generate full-sized matrix from blocks.
MatrixXcd blkdiag(const list<MatrixXcd> & matrix_list, int totalSize) {
    MatrixXcd bdm = MatrixXcd::Zero(totalSize, totalSize);
    int curr_index = 0;
    for (const MatrixXcd & mat : matrix_list)
    {
        if (mat.cols() == 0) {
            continue;
        }
        bdm.block(curr_index, curr_index, mat.rows(), mat.cols()) = mat;
        curr_index += (int) mat.rows();
    }
    return bdm;
}

MatrixXd blkdiag(const list<MatrixXd> & matrix_list, int totalSize) {
    MatrixXd bdm = MatrixXd::Zero(totalSize, totalSize);
    int curr_index = 0;
    for (const MatrixXd & mat : matrix_list)
    {
        if (mat.cols() == 0) {
            continue;
        }
        bdm.block(curr_index, curr_index, mat.rows(), mat.cols()) = mat;
        curr_index += (int) mat.rows();
    }
    return bdm;
}

list<MatrixXcd> blkdiag(const list<list<MatrixXcd>> & matrix_doubleList) {
    list<MatrixXcd> bdmlist;
    for (const list<MatrixXcd> & l : matrix_doubleList) {
        int size = 0;
        for (const MatrixXcd & m : l) {
            size += m.rows();
        }
        bdmlist.emplace_back( blkdiag(l, size) );
    }
    return bdmlist;
}

list<MatrixXd> blkdiag(const list<list<MatrixXd>> & matrix_doubleList) {
    list<MatrixXd> bdmlist;
    for (const list<MatrixXd> & l : matrix_doubleList) {
        int size = 0;
        for (const MatrixXd & m : l) {
            size += m.rows();
        }
        bdmlist.emplace_back( blkdiag(l, size) );
    }
    return bdmlist;
}

// Save vector containing eigenvalues to file.
void saveEnergies(const vector<double> & ergs, const std::string & path) {
    std::ofstream ergFile;
    ergFile.open(path);
    for (double d : ergs) {
        ergFile << d << "\n";
    }
    ergFile.close();
}
