"""
    dpclustgibbs(y::Array{Real, 1}, N::Array{Real, 1}; <keyword arguments>)
Perform dirichlet clustering on the variant allele frequency distribution of cancer sequencing data and find the number of clusters that the data supports, y is a vector of the number of reads reporting each mutant, N is the total depth at each locus.
...
## Arguments
- `iterations = 1000`: number of iterations of the gibbs samples
- `C = 30`: Max number of clusters to consider
- `burninstart = round(Int64, iterations/2)`: Burn in of the gibbs samples
- `bw = 0.01`: Bandwidth of density estimation
- `maxx = 0.7`:
- `cutoffweight = 0.05`: Minimum weight to be called a cluster
- `verbose = true`: Show progress of gibbs sampling with `ProgressMeter` package
...
"""
function dpclustgibbs(y, N;
    iterations = 1000,
    C = 30, #max number of clusters
    burninstart = round(Int64, iterations/2),
    bw = 0.01, # bandwidth of density estimation
    maxx = 0.7,
    cutoffweight = 0.05, #minimum weight to be called a cluster
    verbose = true)

    sum(y .== 0) == 0 || error("Some mutations have VAF = 0.0, make sure these mutations are removed before clustering")

    totalCopyNumber = ones(length(y))
    normalCopyNumber = 2 * ones(length(y))

    # Hyperparameters for alpha, same as Nik-Zainal et al
    A = 0.01
    B = 0.01

    nummuts = length(y)

    # Set up arrays and matrices for recording samples
    π = zeros(iterations, C)
    V = ones(iterations, C)
    S = zeros(Int64, iterations, nummuts)
    PrS = zeros(nummuts, C)
    α = zeros(iterations)
    mutBurdens = zeros(iterations, C, nummuts)

    mutCopyNum = y ./ N

    lower = minimum(mutCopyNum)
    upper = maximum(mutCopyNum)
    difference = upper - lower
    lower = maximum([0.0001, lower - difference/10])
    upper = minimum([upper + difference/10, 0.999])
    # randomise starting positions of clusters
    π[1, :] = rand(Uniform(lower, upper), C)
    for c in 1:C
        mutBurdens[1, c, :] = π[1, c]
    end

    α[1] = 1.0
    V[1, 1:(C - 1)] = 0.5

    if verbose == true
      p = Progress(iterations, 1, "Gibbs sampling progress: ", 30)
    end

    for m in 2:iterations
        @inbounds @simd for k in 1:nummuts
            #Binomial log-likelihood
            PrS[k, 1] = log.(V[m .- 1, 1]) .+ (y[k] .* log.(mutBurdens[m-1, 1, k])) .+
            (N[k] .- y[k]) .* log.(1 .- mutBurdens[m - 1, 1, k])
            PrS[k, 2:C] = allocate(V[m-1, :], mutBurdens[m-1, :, k], y, N, k, 2:C)
            PrS[k, :] = PrS[k, :] .- maximum(PrS[k, :])
            PrS[k, :] = exp.(PrS[k, :])
            PrS[k, :] = PrS[k, :] ./ sum(PrS[k, :])
        end

        S[m, :] = map(k -> sum(rand(Multinomial(1, PrS[k, :])) .* collect(1:length(PrS[k, :]))), 1:nummuts)

        # Update stick-breaking weights
        V[m, 1:(C-1)] = map(h -> rand(Beta(1+sum(S[m, :] .== h), α[m - 1] + sum(S[m, :] .> h))), 1:(C-1))

        V[m, [V[m, 1:(C-1)] .== 1.0; false]] = 0.9999

        countsPerCopyNum = N

        mutBurdens[m, :, :] = mutBurdens[m - 1, :, :]
        @inbounds @simd for c in unique(S[m, :])
          αp = sum(y[S[m, :] .== c])
          βp = 1./sum(countsPerCopyNum[S[m, :] .== c])
          π[m, c] = minimum([rand(Gamma(αp, βp)), 0.999])
          mutBurdens[m, c, :] = π[m, c]
        end

        α[m] = rand(Gamma(C + A - 1, 1/(B - sum(log.(1-V[m, 1:(C-1)])))))

        if verbose == true
          next!(p)
        end

    end

    dp = DPout(S, V, π, α)

    DF, wts = getdensity(dp, iterations; burninstart = burninstart, bw = bw, maxx = maxx)
    wtsout, clonefreq, allwts, allfreq = summariseoutput(dp, wts, iterations; burninstart = burninstart, cutoffweight = cutoffweight)

    sortind = sortperm(clonefreq)
    return DPresults(DF, wts, length(wtsout), wtsout[sortind], clonefreq[sortind], allwts, allfreq, dp, TargetData(y, N, mutCopyNum))
end

function allocate(V, pi, obsy, obsN, currk, jvec)

    out = zeros(length(jvec))

    @inbounds @simd for j in jvec
        out[j-1] = log.(V[j]) .+ sum(log.(1 .- V[1:(j-1)])) .+ obsy[currk] .*log.(pi[j]) .+ (obsN[currk] .- obsy[currk]) .* log.(1-pi[j])
    end

    return out
end

function getdensity(dp, iterations; burninstart = 500, bw = 1.0, maxx = 0.5)

    wts = zeros( size(dp.V)[1], size(dp.V)[2]);
    wts[:, 1] = dp.V[:, 1]
    wts[:, 2] = dp.V[:, 2] .* (1 .- dp.V[:, 1])

    for i in 3:size(wts)[2]
        wts[:, i] = dp.V[:, i] .* prod((1 .- dp.V[:, (1:i .- 1)]), 2)

    end

    postints = zeros(512, iterations -  burninstart + 1)

    xx = kde(dp.π[burninstart - 1, :], weights = wts[burninstart, :]./(sum(wts[burninstart, :])),
                npoints = 512, boundary = (0, maxx), bandwidth = bw).x

    for i in burninstart:iterations
        postints[:, i - burninstart + 1] = kde(dp.π[i - 1, :], weights = wts[i, :]./(sum(wts[i, :])),
        npoints = 512, boundary = (0, maxx), bandwidth = bw).density
    end

    meanv = mean(postints, 2)[:]
    lq = mapslices(x -> quantile(x, 0.025), postints, 2)[:]
    uq = mapslices(x -> quantile(x, 0.975), postints, 2)[:]

    DF = DataFrame(mean = meanv, lq = lq, uq = uq, x = collect(xx))

    return DF, wts
end

function summariseoutput(dp, wts, iterations; burninstart = 1000, cutoffweight = 0.05)

    postwts = wts[burninstart:iterations, :]
    meanwts = mean(postwts, 1)
    clonewts = meanwts[meanwts.>cutoffweight]

    clonefrequency = mean(dp.π[burninstart:iterations, :], 1)
    largeclonefrequency = clonefrequency[meanwts.>cutoffweight]

    clonefrequency = clonefrequency[:]
    meanwts = meanwts[:]

    sortind = sortperm(clonefrequency)

    return clonewts, largeclonefrequency, meanwts[sortind], clonefrequency[sortind]
end
