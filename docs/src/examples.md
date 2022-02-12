# Examples

In this section, we document how the functionality of `RandomMatrixDistributions.jl` can be used to investigate a few standard exampled in random matrix theory.

## Spectrum of a Spiked Wishart 

In this example, we explore the spectral distribution of a spiked Wishart matrix.
This is a matrix of the form ``XX'``, where ``X`` is a ``p \times n`` matrix whose columns are independent ``X \sim \mathcal{N}(0, \Sigma)``, where ``\Sigma = I + H``, where ``H`` is a low-rank symmetric non-negative definite matrix.

The entries of ``X`` can be real or complex, and this is determined by the Dyson parameter ``\beta``.
For this example, we will stick to real entries, setting ``\beta = 1``.

We set our parameters, taking ``H`` to have non-zero eigenvalues of ``0.5, 2, 3``.
These are the "spikes" of the spiked Wishart distribution.

The `SpikedWishart` type encodes the distribution of ``XX'``, and we can compute the corresponding eigenvalue distribution with `EigvalDist`:
Setting `scaled = true` gives the distribution of ``XX'/n``, which is the sample covariance matrix of ``X``.
The scaling is chosen so that the eigenvalue spectrum converges to an appropriate limiting distribution.
```@example specdist
using RandomMatrixDistributions

n = 200
p = 100
beta = 1
spikes = [0.5, 2, 3]

dmat = SpikedWishart(beta, n, p, spikes, scaled = true)
deig = EigvalDist(dmat)
```

We can sample from `deig` as we would from any `Distributions.jl`-style distribution.
In the following, I explicitly supply an object of type `AbstractRNG`, but if left unspecified, the sampler will use `Random.default_rng()`
```@example specdist
using Random
λs = rand(MersenneTwister(0), deig)
```

We can also compute various limiting distributions from `EigvalDist`-type objects.
For example, the `bulk_dist` function computes the limiting bulk (or spectral distribution).

In this case, the limit is taken as ``n, p \rightarrow \infty`` in such a way that ``p/n = \gamma`` is fixed.
```@example specdist
dspec = bulk_dist(deig)
```

To see this behavior more explicitly, let's plot our sampled eigenvalues together with the limiting bulk distribution:
```@example specdist
using Plots, Distributions
histogram(λs, bins = 30, normed = true, legend=false)
plot!(x -> pdf(dspec, x))
```

The two eigenvalues outside the bulk are the supercritical spikes.
Since the smallest spike of ``0.5`` is smaller than ``\sqrt{\gamma}``, it doesn't emerge from the bulk and is called subcritical.

The `supercrit_dist` function can be used to compute the approximate distribution of the supercritical eigenvalues.
In this case, they are independent Gaussians, which are biased versions of the corresponding population spikes of ``2`` and ``3``:
```@example specdist
supercrit_dist(deig)
```


## Largest eigenvalues of a GUE

The largest eigenvalue of a GUE matrix has a Tracy-Widom limiting distribution.
In particular, if ``\lambda_1`` is the largest eigenvalue of an ``n \times n`` GUE matrix, then ``n^{2/3}(\lambda_1 - 2) \stackrel{\mathrm{d}}{\rightarrow} \mathrm{TW}(2)``, where the ``2`` in ``\mathrm{TW}(2)`` denotes the Dyson parameter, as is appropriate for the complex entries of the GUE.

To repeatedly sample ``\lambda_1`` values, we use the `randeigstat` function.
It takes as arguments a random matrix distribution (a GUE matrix is a Wigner matrix with ``\beta = 2`` and no spikes), a function of a matrix argument, and the number of samples.

For the function to work correctly, the function of a matrix argument must be an eigenvalue statistic - that it, it should  depend only on the eigenvalues of that matrix.
For example, `eigmax` computes the largest eigenvalue of a matrix, and so is an eigenvalue statistic.

This is implemented in this way rather than explicitly requiring that the eigenvalue statistic be specified explicitly as a function of the eigenvalues, since often this would lead to unnecessary computation.
For example, the trace is an eigenvalue statistic, but can be evaluated very quickly without computing eigenvalues.

In the following, we sample 1000 realizations of ``n^{2/3}(\lambda_1 - 2)``:
```@example max_eig
using RandomMatrixDistributions

n = 500
reps = 1000

d = SpikedWigner(2, n, scaled = true)
λs = randeigstat(d, eigmax, reps)

λs_scaled = n^(2/3) * (λs .- 2)
```

`RandomMatrixDistributions` implements `TracyWidom(β)` according to the `Distributions.jl`-style API.
We can plot its pdf to compare to the histogram of the scaled ``\lambda_1``s using the `pdf` function inherited from that package.
```@example max_eig
using Plots, Distributions

histogram(λs_scaled, normed = true, bins = 100)
plot!(x -> pdf(TracyWidom(2), x))
```

Other relevant methods from `Distributions.jl` are also implemented, including the `quantile` function.
In the following, we make use of this in order to produce a QQ-plot of our sample against the Tracy-Widom distribution.

To this end, we evaluate the Tracy-Widom quantiles at `k/(reps + 1)` for `k = 1, ..., reps`:
```@example max_eig
quantile_positions = (1:reps)/(reps+1)
TW_quantiles = quantile(TracyWidom(2), quantile_positions)
nothing # hide
```

The QQ-plot is then just a plot of our sorted sample against these theoretical quantile.
Adding a red line with slope 1 and intercept 0 demonstrates the close agreements between the theoretical and sample quantiles:
```@example max_eig
scatter(
    TW_quantiles, sort(λs_scaled),
    xlabel = "TW(2) theoretical quantiles",
    ylabel = "Scaled largest eigenvalues",
    legend=false, markersize=2, seriescolor=nothing)

plot!(x->x)
```

## Speed comparisons for eigenvalue computations

```@setup timing
using RandomMatrixDistributions, LinearAlgebra
dist = SpikedWigner(1, 500)
@time eigvals.(rand(dist))
@time rand(EigvalDist(dist))
dist = SpikedWigner(1, 500, [1, 1, 1, 1])
@time eigvals.(rand(dist))
@time rand(EigvalDist(dist))
```

Sampling from `EigvalDist`-type distributions is done internally by sampling a banded matrix whose eigenvalues have the same distribution as those of a corresponding dense matrix, and using the functionality of `BandedMatrices.jl` to efficiently compute the eigenvalues of the resulting matrix.

In te following example, we compare the resources used for sampling a dense ``500 \times 500`` GOE matrix, and then computing its eigenvalues to those used when directly sampling from the corresponding `EigvalDist`-type distribution:
```@example timing
using RandomMatrixDistributions, LinearAlgebra
dist = SpikedWigner(1, 500)
@time eigvals.(rand(dist, 200))
@time rand(EigvalDist(dist), 200)
nothing # hide
```

Julia computes the eigenvalues of Hermitian matrices efficiently by reducing to a tridiagonal form, so there sampling a tridiagonal matrix directly improves the running time by a constant factor.
On the other hand, there is a serious gain in terms of allocations, as the representation of a banded matrix requires ``O(n)`` in comparison to the ``O(n^2)`` required for a dense matrix.

We demonstrate this by sampling eigenvalues of GOE matrices using both methods with increasing values of ``n``.
The results are plotted on a log-log plot, since resourse use in all cases is of the form `A n^\kappa`:
```@example timing
nvals = 50:50:1000
reps = 100

times = Array{Float64}(undef, length(nvals), 2)
allocs = Array{Float64}(undef, length(nvals), 2)

for (i, n) in enumerate(nvals)
    d = SpikedWigner(1, n)
    
    dense = @timed eigvals.(rand(d, reps))
    times[i,1] = dense.time
    allocs[i,1] = dense.bytes
    
    banded = @timed rand(EigvalDist(d), reps)
    times[i,2] = banded.time
    allocs[i,2] = banded.bytes
end
```

```@example timing
using Plots # hide
p1 = plot(nvals, times, # hide
    label = ["Dense matrix" "Banded matrix"], # hide
    xlabel = "n", ylabel = "Elapsed time (seconds)", # hide
    legend = :bottomright, # hide
    xaxis=:log, yaxis=:log # hide
) # hide

p2 = plot(nvals, allocs, # hide
    label = ["Dense matrix" "Banded matrix"], # hide
    xlabel = "n", ylabel = "Allocations (bytes)", # hide
    legend = false, # hide
    xaxis=:log, yaxis=:log # hide
) # hide

plot(p1, p2) # hide
```

For matrices with ``r`` spikes, the running time and number of allocations required to compute the eigenvalues of a banded representation scales linearly in ``r``, as we see in the following:
```@example timing
spikenums = 2:2:50
n = 100
reps = 200

times = Array{Float64}(undef, length(spikenums), 2)
allocs = Array{Float64}(undef, length(spikenums), 2)

for (i, numspikes) in enumerate(spikenums)
    d = SpikedWigner(1, n, ones(numspikes), scaled=false)
    
    dense = @timed eigvals.(rand(d, reps))
    times[i,1] = dense.time
    allocs[i,1] = dense.bytes
    
    banded = @timed rand(EigvalDist(d), reps)
    times[i,2] = banded.time
    allocs[i,2] = banded.bytes
end
```

```@example timing
p1 = plot(spikenums, times, # hide
    label = ["Dense matrix" "Banded matrix"], # hide
    xlabel = "Number of spikes", ylabel = "Elapsed time (seconds)", # hide
    legend=false, # hide
    xaxis=:log, yaxis=:log # hide
) # hide

p2 = plot(spikenums, allocs, # hide
    label = ["Dense matrix" "Banded matrix"], # hide
    xlabel = "Number of spikes", ylabel = "Allocations (bytes)", # hide
    legend = :bottomright, # hide
    xaxis=:log, yaxis=:log # hide
) # hide

plot(p1, p2) # hide
```
