RandomMatrixDistributions.jl
============================
[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://damian-t-p.github.io/RandomMatrixDistributions.jl/dev/)
[![Build Status](https://github.com/damian-t-p/RandomMatrixDistributions.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/damian-t-p/RandomMatrixDistributions.jl/actions?query=workflow)
[![codecov](https://codecov.io/gh/damian-t-p/RandomMatrixDistributions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/damian-t-p/RandomMatrixDistributions.jl)

A Julia package containing `Distributions.jl`-type specifications for various distributions arising from random matrix theory.

# Currently implemented distributions

## Matrix distributions

* `SpikedWigner(beta, n, spikes; scaled=false)`: Wigner distribution with an added spike matrix.

	`spikes` is an array `[s1, ..., sr]` such that the diagonal matrix with diagonal sqrt(n)(*s*<sub>1</sub>, ... , *s*<sub>r</sub>, 0, ..., 0) is added to a white Wigner matrix.

* `SpikedWishart(beta, n, p, spikes; scaled=false)`: Wishart distribution with spiked covariance [1].

	`spikes` is an array `[s1, ..., sr]` such that the Wishart covariance is diagonal with entries  (1 + *s*<sub>1</sub>, ... , 1 + *s*<sub>r</sub>, 1, ..., 1).
	
* `Jacobi(beta, n1, n2, p)`: Random matrices of the form *E*(*E*+*H*)<sup>-1</sup>. Here *E* and *H* are (*n*<sub>1</sub>, *p*) and (*n*<sub>2</sub>, *p*) white Wisharts respectively. [2]

Specifying `scaled=true` in `SpikedWigner` and `SpikedWishart` scales the matrices by an appropriate function of *n* so that the corresponding bulks converge to the semicircle and Marchenko-Pastur laws respectively.
Due to the inverse in the definition of the Jacobi ensemble, no scaling is necessary for `Jacobi`,

Normal entries in Gaussian ensembles are scaled to have variance 1.

## Limiting eigenvalue distributions

The package implements the following types:
* `MarchenkoPastur(gamma)`: Limiting empirical spectral density of a real white Wishart matrix with *p*/*n* -> *gamma* as long as 0 < *gamma* < 1.
* `TracyWidom(beta)`: Limiting distribution of the maximum eigenvalue of many random matrix ensembles with Dyson parameter beta [3].
* `Wachter(gamma1, gamma2)`: Limiting empirical spectral density of *S*<sub>1</sub> *S*<sub>2</sub><sup>-1</sup>. Here *S*<sub>1</sub> and *S*<sub>2</sub> are sample covariance matrices with *n*<sub>1</sub>/*p* -> *gamma*<sub>1</sub> and *n*<sub>2</sub>/*p* -> *gamma*<sub>2</sub>.

It also implements the following functions for computing eigenvalue distributions:
* `supercrit_dist(E)`: Approximate distribution of the supercritical eigenvalues of a matrix drawn from the ensemble E.
Currently implemented for Wishart with beta = 1 [4], beta = 2 [5] and Wigner with beta = 1 [6], beta = 2 [7]

# Efficient samplers
   The function `randeigvals` efficiently samples from the distribution of eigenvalues of the implemented random matrix distributions. It does this by generating a tridiagonal or banded matrix with eigenvalue equal in distribution to the specified model. 

# Examples
   See the documentation.

# References
[1] Dumitriu & Edelman, "Matrix Models for beta ensembles," Journal of Mathematical Physics, 11 (2002).

[2] Killip & Nenciu, "Matrix Models for Circular Ensembles," International Mathematics Research Notices, 50 (2004).

[3] Bornemann, "On the numerical evaluation of distributions in random matrix theory: a review," (2010).

[4] Baik, Ben Arous & Peche, "Phase transition of the largest eigenvalue for nonnull complex sample covariance matrices," The Annals of Probability, 33 (2005).

[5] Paul, "Asymptotics of sample eigenstructure for a large dimensional spiked covariance model," Statistica Sinica, 17 (2007).

[6] Feral, Peche, "The largest eigenvalue of rank one deformation of large wigner matrices," Commun. Math. Phys., 272 (2007).

[7] Peche, "The largest eigenvalue of small rank perturbations of Hermitian random matrices," Probability Theory and Related Fields, 134 (2006).

