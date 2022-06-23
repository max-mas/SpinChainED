# SpinChainED
This project models a zig-zag triangular spin ladder or, equivalently, a Heisenberg Chain with next-nearest-neighbour coupling.

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

Required: Eigen, OpenMP, Intel MKL
