var documenterSearchIndex = {"docs":
[{"location":"library/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Pages = [\"library.md\"]","category":"page"},{"location":"library/#Distributions","page":"Library","title":"Distributions","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"EigvalDist\nJacobi\nSpikedWigner\nSpikedWishart","category":"page"},{"location":"library/#RandomMatrixDistributions.EigvalDist","page":"Library","title":"RandomMatrixDistributions.EigvalDist","text":"EigvalDist(matrixdist::MatrixDistribution)\n\nThe joint distribution of eigenvalues sampled from a random matrix with the distribution matrixdist.\n\n\n\n\n\n","category":"type"},{"location":"library/#RandomMatrixDistributions.Jacobi","page":"Library","title":"RandomMatrixDistributions.Jacobi","text":"Jacobi(β::Int, n₁::Int, n₂::Int, p::Int)\n\nDistribution of a p×p Jacobi matrix.\n\nIf E ~ Wishartₚ(I, n₁) and H ~ Wishartₚ(I, n₂) are independent with Dyson parameter β, then E(E + H)⁻¹ has Jacobi(β, n₁, n₂, p) distribution.\n\nIf λᵢ are the eigenvalues of EH⁻¹ and μᵢ are the Jacobi eigenvalues, then μᵢ = λᵢ/(1 + λᵢ) and λᵢ = μᵢ/(1 - μᵢ).\n\n\n\n\n\n","category":"type"},{"location":"library/#RandomMatrixDistributions.SpikedWigner","page":"Library","title":"RandomMatrixDistributions.SpikedWigner","text":"SpikedWigner(β::Int, n::Int[, spikes::Vector{Float64}, scaled::Bool=false])\n\nDistribution on an n×n spiked Gaussian Wigner matrix.\n\nWigner matrices are Hermitian with independent real, complex or quaternion standard Gaussian entries depending on whether β = 1, 2 or 4.\n\nA diagonal matrix with entries given by spikes multiplied by √n is added to produce the spiked Wigner matrix.\n\nIf scaled == true, then the resulting matrix is divided by √n so that its bulk distribution converges to the semicircle law supported on [-2, 2].\n\n\n\n\n\n","category":"type"},{"location":"library/#RandomMatrixDistributions.SpikedWishart","page":"Library","title":"RandomMatrixDistributions.SpikedWishart","text":"SpikedWishart(β::Int, n::Int, p::Int[, spikes::Vector{Float64}, scaled::Bool=false)\n\nDistribution of a p×p spiked Wishart matrix.\n\nIf X is a p×n matrix with independent real, complex or quaternion standard Gaussian entries, depending on whether β = 1, 2 or 4, then XX† has a Wishart(β, n, p) distribution.\n\nIf Λ is a diagonal matrix whose entries are √(1 .+ spikes), then ΛXX†Λ has a SpikedWishart(β, n, p, spikes) distribution.\n\nIf scaled == true, then the resulting matrix is divided by p so that its bulk distribution converges to the Marchenko-Pastur law.\n\n\n\n\n\n","category":"type"},{"location":"library/#Methods-for-random-matrix-ensembles","page":"Library","title":"Methods for random matrix ensembles","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"bulk_dist\nrandeigvals\nrandeigstat\nsupercrit_dist","category":"page"},{"location":"library/#RandomMatrixDistributions.bulk_dist","page":"Library","title":"RandomMatrixDistributions.bulk_dist","text":"bulk_dist(d::Union{MatrixDistribution, EigvalDist})\n\nCompute the limiting spectral distribution of d.\n\n\n\n\n\n","category":"function"},{"location":"library/#RandomMatrixDistributions.randeigvals","page":"Library","title":"RandomMatrixDistributions.randeigvals","text":"randeigvals([rng::AbstractRNG, ]d::MatrixDistribution)\n\nSample a vector of eigenvalues of a matrix drawn from the matrix ensemble d.\n\n\n\n\n\n","category":"function"},{"location":"library/#RandomMatrixDistributions.randeigstat","page":"Library","title":"RandomMatrixDistributions.randeigstat","text":"randeigstat([rng::AbstractRNG, ]d::MatrixDistribution, eigstat::Function, n::Int)\n\nSample n realisations of the eigenvalue statistic eigstat evaluated at a matrices drawn from the ensemble d.\n\neigstat is a function of a square matrix argument whose value depends only on the eigenvalues of that matrix.\n\nUsage\n\nrangeigstat(SpikedWigner(2, 50), eigmax, 100)\n\n\n\n\n\n","category":"function"},{"location":"library/#RandomMatrixDistributions.supercrit_dist","page":"Library","title":"RandomMatrixDistributions.supercrit_dist","text":"supercrit_dist(d::Union{MatrixDistribution, EigvalDist})\n\nCompute the approximate joint distribution of the supercritical eigenvalues of the ensemble d.\n\n\n\n\n\n","category":"function"},{"location":"library/#Limiting-eigenvalue-densities","page":"Library","title":"Limiting eigenvalue densities","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"MarchenkoPastur\nTracyWidom\nWachter","category":"page"},{"location":"library/#RandomMatrixDistributions.MarchenkoPastur","page":"Library","title":"RandomMatrixDistributions.MarchenkoPastur","text":"MarchenkoPastur(γ::Real)\n\nMarchenko-Pastur distribution, where 0 < γ ≤ 1.\n\nThe limiting spectral distribution of a p×p covariance matrix of n standard normal observations, where p/n → γ.\n\n\n\n\n\n","category":"type"},{"location":"library/#RandomMatrixDistributions.TracyWidom","page":"Library","title":"RandomMatrixDistributions.TracyWidom","text":"TracyWidom(β::Int)\n\nTracy-Widom distribution with Dyson parameter β.\n\nThe limiting distribution of the largest eigenvalue of a GOE (β = 1), GUE (β = 2) or GSE (β = 4) matrix.\n\n\n\n\n\n","category":"type"},{"location":"library/#RandomMatrixDistributions.Wachter","page":"Library","title":"RandomMatrixDistributions.Wachter","text":"Wachter(γ₁::Real, γ₂::Real)\n\nWachter distribution, where 0 ≤ γ₂ < 1.\n\nLet Σ₁ and Σ₂ be p×p covariance matrices of n₁ and n₂ standard normal observations respectively.\n\nIf p/n₁ → γ₁ and p/n₂ → γ₂, then Σ₁ Σ₂⁻¹ has a limiting spectral distribution of Wachter(γ₁, γ₂).\n\n\n\n\n\n","category":"type"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Spectrum-of-a-Spiked-Wishart","page":"Examples","title":"Spectrum of a Spiked Wishart","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this example, we explore the spectral distribution of a spiked Wishart matrix. This is a matrix of the form XX, where X is a p times n matrix whose columns are independent X sim mathcalN(0 Sigma), where Sigma = I + H, where H is a low-rank symmetric non-negative definite matrix.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The entries of X can be real or complex, and this is determined by the Dyson parameter beta. For this example, we will stick to real entries, setting beta = 1.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We set our parameters, taking H to have non-zero eigenvalues of 05 2 3. These are the \"spikes\" of the spiked Wishart distribution.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The SpikedWishart type encodes the distribution of XX, and we can compute the corresponding eigenvalue distribution with EigvalDist: Setting scaled = true gives the distribution of XXn, which is the sample covariance matrix of X. The scaling is chosen so that the eigenvalue spectrum converges to an appropriate limiting distribution.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using RandomMatrixDistributions\n\nn = 200\np = 100\nbeta = 1\nspikes = [0.5, 2, 3]\n\ndmat = SpikedWishart(beta, n, p, spikes, scaled = true)\ndeig = EigvalDist(dmat)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We can sample from deig as we would from any Distributions.jl-style distribution. In the following, I explicitly supply an object of type AbstractRNG, but if left unspecified, the sampler will use Random.default_rng()","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Random\nλs = rand(MersenneTwister(0), deig)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We can also compute various limiting distributions from EigvalDist-type objects. For example, the bulk_dist function computes the limiting bulk (or spectral distribution).","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this case, the limit is taken as n p rightarrow infty in such a way that pn = gamma is fixed.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"dspec = bulk_dist(deig)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To see this behavior more explicitly, let's plot our sampled eigenvalues together with the limiting bulk distribution:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Plots, Distributions\nhistogram(λs, bins = 30, normed = true, legend=false)\nplot!(x -> pdf(dspec, x))","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The two eigenvalues outside the bulk are the supercritical spikes. Since the smallest spike of 05 is smaller than sqrtgamma, it doesn't emerge from the bulk and is called subcritical.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The supercrit_dist function can be used to compute the approximate distribution of the supercritical eigenvalues. In this case, they are independent Gaussians, which are biased versions of the corresponding population spikes of 2 and 3:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"supercrit_dist(deig)","category":"page"},{"location":"examples/#Largest-eigenvalues-of-a-GUE","page":"Examples","title":"Largest eigenvalues of a GUE","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The largest eigenvalue of a GUE matrix has a Tracy-Widom limiting distribution. In particular, if lambda_1 is the largest eigenvalue of an n times n GUE matrix, then n^23(lambda_1 - 2) stackrelmathrmdrightarrow mathrmTW(2), where the 2 in mathrmTW(2) denotes the Dyson parameter, as is appropriate for the complex entries of the GUE.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To repeatedly sample lamdda_1 values, we use the randeigstat function. It takes as arguments a random matrix distribution (a GUE matrix is a Wigner matrix with beta = 2 and no spikes), a function of a matrix argument, and the number of samples.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"For the function to work correctly, the function of a matrix argument must be an eigenvalue statistic - that it, it should  depend only on the eigenvalues of that matrix. For example, eigmax computes the largest eigenvalue of a matrix, and so is an eigenvalue statistic.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This is implemented in this way rather than explicitly requiring that the eigenvalue statistic be specified explicitly as a function of the eigenvalues, since often this would lead to unnecessary computation. For example, the trace is an eigenvalue statistic, but can be evaluated very quickly without computing eigenvalues.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In the following, we sample 1000 realizations of n^23(lambda_1 - 2):","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using RandomMatrixDistributions\n\nn = 500\nreps = 1000\n\nd = SpikedWigner(2, n, scaled = true)\nλs = randeigstat(d, eigmax, reps)\n\nλs_scaled = n^(2/3) * (λs .- 2)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"RandomMatrixDistributions implements TracyWidom(β) according to the Distributions.jl-style API. We can plot its pdf to compare to the histogram of the scaled lambda_1s using the pdf function inherited from that package.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Plots, Distributions\n\nhistogram(λs_scaled, normed = true, bins = 100)\nplot!(x -> pdf(TracyWidom(2), x))","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Other relevant methods from Distributions.jl are also implemented, including the quantile function. In the following, we make use of this in order to produce a QQ-plot of our sample against the Tracy-Widom distribution.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To this end, we evaluate the Tracy-Widom quantiles at k/(reps + 1) for k = 1, ..., reps:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"quantile_positions = (1:reps)/(reps+1)\nTW_quantiles = quantile(TracyWidom(2), quantile_positions)\nnothing # hide","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The QQ-plot is then just a plot of our sorted sample against these theoretical quantile. Adding a red line with slope 1 and intercept 0 demonstrates the close agreements between the theoretical and sample quantiles:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"scatter(\n    TW_quantiles, sort(λs_scaled),\n    xlabel = \"TW(2) theoretical quantiles\",\n    ylabel = \"Scaled largest eigenvalues\",\n    legend=false, markersize=2, seriescolor=nothing)\n\nplot!(x->x)","category":"page"},{"location":"#RandomMatrixDistributions.jl","page":"Home","title":"RandomMatrixDistributions.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package containing Distributions.jl-type specifications for various distributions arising in random matrix theory with a focus on spiked ensembles.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"RandomMatrixDistributions is a registered package and can be installed by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pkg.add(\"RandomMatrixDistributions\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"in the Julia REPL.","category":"page"}]
}
