using DPclustering
using Base.Test

tests = ["nikzainal"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
