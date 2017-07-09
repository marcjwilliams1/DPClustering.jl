type DPout
    S
    V
    π
    α
end

type DPresults
    DF::DataFrame
    weights
    nclones
    cloneweights
    clonefrequencies
    posterior::DPout
end
