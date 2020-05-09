RandomMatrixDistributions.jl
============================

[![Build Status](https://travis-ci.com/damian-t-p/RandomMatrixDistributions.jl.svg?branch=master)](https://travis-ci.com/damian-t-p/RandomMatrixDistributions.jl)
[![codecov](https://codecov.io/gh/damian-t-p/RandomMatrixDistributions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/damian-t-p/RandomMatrixDistributions.jl)

A Julia package containing =Distribution.jl=-type specifications for various distributions arising from random matrix theory.

# Currently implemented distributions

## Matrix distributions

* `SpikedWishart(beta, n, p, spikes; scaled=false)`: Wishart distribution with spiked covariance (sampler for more than one spike implemented only for the real case. [1]

	`spikes` is an array `[s1, ..., sr]` such that the Wishart covariance is diagonal with entries  *s<sub>1</sub>, ... , s<sub>r</sub>, 1, \dotsc, 1*.
	
* `Jacobi(beta, n1, n2, p)`: Random matrices of the form $E(E+H)^{-1}$. Here $E$ and *H* are *(n<sub>1</sub>, p)* and *(n<sub>2</sub>, p)* white Wisharts respectively. [2]

## Limiting eigenvalue distributions
* `MarchenkoPastur(gamma)`: Limiting empirical spectral density of a real white Wishart matrix with *p/ -> \gamma/* as long as *0 < &gamma < 1*.
* `TracyWidom(beta)`: Limiting distribution of the maximum eigenvalue of many random matrix ensembles with Dyson parameter beta.
* `Wachter(gamma1, gamma2)`: Limiting empirical spectral density of $\Sigma_{1} \Sigma_{2}^{-1}$. Here $\Sigma_{1}$ and $\Sigma_{2}$ are sample covariance matrices with $n_1/p \rightarrow \gamma_1$ and $n_{2}/p \rightarrow \gamma_{2}$.

# Efficient samplers
   The function `randeigvals` efficiently samples from the distribution of eigenvalues of the implemented random matrix distributions. It does this by generating a tridiagonal or banded matrix with eigenvalue equal in distribution to the specified model. 

# Examples
   An Jupyter notebook demonstrating all of the implemented eigenvalue samplers is provided in `/examples/eigenvalue-simulation.ipynb`.

# References
[1] Dumitriu & Edelman, Matrix Models for beta ensembles, Journal of Mathematical Physics, (11), (2002).

[2] Killip & Nenciu, Matrix Models for Circular Ensembles, International Mathematics Research Notices, 50, (2004).

