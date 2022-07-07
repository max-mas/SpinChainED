#include "main.h"

using std::vector;

int main(int argc, char* argv[]) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    const int cpu_cnt = (int) std::thread::hardware_concurrency() / 2;
    omp_set_num_threads(cpu_cnt);

    int minN = 18;
    int maxN = 18;
    std::cout << "DQT Specific Heat, 5000 Data Points:" << std::endl;
    for (int N = minN; N <= maxN; N+=2) {
        const Eigen::SparseMatrix<std::complex<double>> S2 = spinOp2_momentum_sparse(N);
        saveSusceptibilityForVaryingTemp_DQT_avg(N, 5000, 0, 50, S2, "", 1);

    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << "[ms]" << std::endl;
    return 0;
}
