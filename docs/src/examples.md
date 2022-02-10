# Examples

## Empirical spectral distributions

In this example, we explore the spectral distribution of a spiked Wishart matrix.
This is a matrix of the form ``XX'``, where ``X`` is a ``p \times n`` matrix whose columns are independent ``X \sim \mathcal{N}(0, \Sigma)``, where ``\Sigma = I + H``, where ``H`` is a low-rank symmetric non-negative definite matrix.

The entries of ``X`` can be real or complex, and this is determined by the Dyson parameter ``\beta``.
For this example, we will stick to real entries, setting ``\beta = 1``.

We set our parameters, taking ``H`` to have non-zero eigenvalues of ``0.5, 2, 3``.
These are the "spikes" of the spiked Wishart distribution.

The `SpikedWishart` type encodes the distribution of ``XX'``, and we can compute the corresponding eigenvalue distribution with ``EigvalDist``:
Setting `scaled = true` gives the distribution of ``XX'/n``, which is the sample covariance matrix of ``X``.
The scaling is chosen so that the eigenvalue spectrum converges to an appropriate limiting distribution.
```jldoctest specdist
using RandomMatrixDistributions

n = 200
p = 100
beta = 1
spikes = [0.5, 2, 3]

dmat = SpikedWishart(beta, n, p, spikes, scaled = true)
deig = EigvalDist(dmat)

# output

EigvalDist(
matrixdist: SpikedWishart(beta=1, n=200, p=100, spikes=[0.5, 2.0, 3.0], scaled=true)
)
```

We can sample from `deig` as we would from any `Distributions.jl`-style distribution.
In the following, I explicitly supply an object of type `AbstractRNG`, but if left unspecified, the sampler will use `Random.default_rng()`
```@meta specdist
DocTestFilters = [
r"0\.1869251107737921(?s).*2\.4313768183766524",
r"[0-9]{5}(\n|$)"
]
```

```jldoctest specdist
using Random
λs = rand(MersenneTwister(0), deig)

# output

100-element Vector{Float64}:
 0.09955938791850548
 0.10653785238288815
 0.12269801419629447
 0.1330483719251576
 0.1392930357410377
 0.15464459410357262
 0.17043480014008375
 0.1869251107737921
 ⋮
 2.4313768183766524
 2.5011406257279383
 2.6929816190365417
 2.762519332554983
 2.870429730291467
 3.914014518969469
 4.507428322596304
```

We can also compute various limiting distributions from `EigvalDist`-type objects.
For example, the `bulk_dist` function computes the limiting bulk (or spectral distribution).

In this case, the limit is taken as ``n, p \rightarrow \infty`` in such a way that ``p/n = \gamma`` is fixed.
```jldoctest specdist
dspec = bulk_dist(deig)

# output

MarchenkoPastur(gamma=0.5)
```

To see this behavior more explicitly, let's plot our sampled eigenvalues together with the limiting bulk distribution:
```@setup specdist_eg
using RandomMatrixDistributions, Random
n = 200
p = 100
beta = 1
spikes = [0.5, 2, 3]
dmat = SpikedWishart(beta, n, p, spikes, scaled = true)
deig = EigvalDist(dmat)
λs = rand(MersenneTwister(0), deig)
dspec = bulk_dist(deig)
```
```@example specdist_eg
using Plots, Distributions
histogram(λs, bins = 30, normed = true, legend=false)
plot!(x -> pdf(dspec, x))
```

The two eigenvalues outside the bulk are the supercritical spikes.
Since the smallest spike of ``0.5`` is smaller than ``\sqrt{\gamma}``, it doesn't emerge from the bulk and is called subcritical.

The `supercrit_dist` function can be used to compute the approximate distribution of the supercritical eigenvalues.
In this case, they are independent Gaussians, which are biased versions of the corresponding population spikes of ``2`` and ``3``:
```jldoctest specdist
supercrit_dist(deig)

# output

DiagNormal(
dim: 2
μ: [3.75, 4.666666666666667]
Σ: [0.039375000000000014 0.0; 0.0 0.07555555555555556]
)
```
