# SpinChainED
This project models a zig-zag triangular spin ladder or, equivalently, a Heisenberg Chain with next-nearest-neighbour coupling, which has the Hamiltonian

![equation](https://latex.codecogs.com/png.image?\large&space;\dpi{120}\color{white}H&space;=&space;J_2\left(\sum_{i=1}^N\vec{S}_i\cdot\vec{S}_{i&plus;1}&space;&plus;&space;\frac{J_1}{J_2}\sum_{i=1}^N\vec{S}_i\cdot\vec{S}_{i&plus;2}&space;\right))

The full spectrum can be computed via ED for given coupling constants and a set system size using a variety of symmetries.
These are:

	- No symmetries, aka naive
	- Magnetization or total spin conservation
	- Translational lattice symmetry using momentum states
	- Parity symmetry
	- Spin-inversion symmetry for the m = S_z = 0 block
	
The spectrum can then be used to calculate specific heat and magnetic susceptibility curves.

These quantities can also be computed using a Dynamical Quantum Typicality (DQT) approach, where thermal expectation values are approximated using scalar products of operators and random vectors. The temperature dependece is generated using an imaginary-time Schrödinger equation iterated via the common RK4 scheme.

The python files in the root of the project can be used to generate plots. The file spinGapFit also provides functionality for approximating the spin energy gap and excitation energy of the system using the thermodynamic quantities computed using DQT.

When building the project using cmake, make sure to adjust the paths of eigen and MKL in CMakeLists.txt. Also make sure to adjust the paths where data files will be saved to a proper location on disk in main.cpp an the plotting files.

When running the program, make sure to use the correct command line arguments. These will be updated to contain DQT functionality in the future.
At the moment, these are:

	- 0: nMin, minimum n for which computations will be performed
	- 1: nMax, minimum n for which computations will be performed
	- 2: dataPointNum, number of data points used in desired computation
	- 3: J_ratios_path, path to a file named J_ratios.txt that contains the values of J_1/J_2 for which computations will be performed
	- 4: Ts, path to a file named Ts.txt that contains the values of T or beta for which computations will be performed
	- 5: isBeta, 0 or 1, indicating whether computations should be carried out in T or beta space
	- 6: start, value of T, beta or J_1/J_2 at which computations will begin
	- 7: end, value of T, beta or J_1/J_2 at which computations will stop
	- 8: flags, seven-digit series of 0s and 1s indicating whether or not the following quantities should be included in computations:
		- spin gap or excitation energies (edit main.cpp to set)
		- ground stat energy per spin
		- specific heat
		- specific heat for varying J_1/J_2
		- susceptibility
		- susceptibility for varying J_1/J_2
		- k-dispersion of states
	- 9: save_to_path, path where data will be saved, must either have the correct directory structure or main.cpp must be edited


Required for build: Eigen, OpenMP (usually part of gcc), Intel MKL

The Hamiltonian-generating methods of this project have been inspired by Anders W. Sandvik, Computational Studies of Quantum Spin Systems,  https://doi.org/10.1063/1.3518900
