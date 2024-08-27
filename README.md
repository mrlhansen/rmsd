# Spectral densities
This code can be used for calculating smeared spectral densities from lattice QCD correlators. This is done by solving an inverse Laplace problem, as described in the paper:

> Martin Hansen, Alessandro Lupo, Nazario Tantalo, On the extraction of spectral densities from lattice correlators, Phys.Rev. D99 (2019) 094508, [arXiv:1903.06476](https://arxiv.org/abs/1903.06476).

To properly understand what the code does, it is advised to carefully read this paper. Moreover, if you decide to use this code, or the implemented algorithm, please cite the paper.

## Usage
Compiling the code requires a working version of GCC and two multi-precision libraries. These libraries are needed because the algorithm solves a linear system with an extremely ill-conditioned matrix.

* [MPFR](https://www.mpfr.org)
* [FLINT](https://flintlib.org/)

When compiling the code, two main programs are produced.

The first program *fakecorr* generates a synthetic correlator from the model described in the previously mentioned paper, and uses the algorithm to reconstruct the spectral density. Because the benchmark model is known exactly, in this case it is possible to compare the exact spectral density to the one reconstructed by the algorithm. When running the program with the option `-nms 1` no noise is added to the correlator.

The second program *realcorr* reads a correlator from a file, and then performs the reconstruction of the spectral density. The input file is expected to be a text file where each row contains a correlator C(t), with each element being delimited by a space, and the different rows are regarded as independent measurements of the correlator. The time extent and the number of measurements must be correctly specified via the command line options. The code assumes that the input correlator is symmetric (i.e. with periodic boundary conditions in time) and that the number of measurements is larger than one.

When running the programs without any arguments, they will show a list of all possible command line options. The different choices for the smearing kernel are defined in `docs/kernels.pdf`.

## Output
The code writes several output files in a separate folder, whose name is the current unix timestamp. The format of the different output files are described in the file `docs/output.pdf`.
