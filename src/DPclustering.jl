module DPClustering

# packages
using DataFrames
using ProgressMeter
using KernelDensity
using Plots
using Distributions
using Colors

import Base.show

export dpclustgibbs,
DPresults,
plotresults,
show

include("types.jl")
include("plotting.jl")
include("util.jl")
include("clustering.jl")

end # module
