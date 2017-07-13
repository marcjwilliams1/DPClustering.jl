type DPout
    S
    V
    π
    α
end

type TargetData

    y::Array{Int64, 1}
    N::Array{Int64, 1}
    VAF::Array{Float64, 1}
    DF::DataFrame

    TargetData(y, N, VAF) =
    new(y, N, VAF, DataFrame(y = y, N = N, VAF = VAF))
end

type DPresults
    DF::DataFrame
    weights
    nclones
    cloneweights
    clonefrequencies
    allcloneweights
    allclonefrequencies
    posterior::DPout
    data::TargetData
end
