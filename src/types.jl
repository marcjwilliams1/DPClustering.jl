mutable struct DPout
    S
    V
    π
    α
end

mutable struct TargetData

    y::Array{Int64, 1}
    N::Array{Int64, 1}
    VAF::Array{Float64, 1}
    DF::DataFrame

    TargetData(y, N, VAF) =
    new(y, N, VAF, DataFrame(y = y, N = N, VAF = VAF))
end

mutable struct DPresults
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
