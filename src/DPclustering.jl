module DPclustering

# packagesusing Distributions
using RCall
using DataFrames
using ProgressMeter
using KernelDensity
using CancerSeqSim
using Gadfly

export dpclustgibbs,
DPresults,
plotresults

include("types.jl")
include("clustering.jl")
include("plotting.jl")
include("ggplotting.jl")

end # module
