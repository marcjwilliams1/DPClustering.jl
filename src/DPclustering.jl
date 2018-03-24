module DPclustering

# packages
using DataFrames
using ProgressMeter
using KernelDensity
using Gadfly
using Distributions
using Colors

import Base.show

export dpclustgibbs,
DPresults,
plotresults,
show

include("types.jl")
include("clustering.jl")
include("plotting.jl")
include("util.jl")

end # module
