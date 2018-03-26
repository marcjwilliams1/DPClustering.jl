#Running DP clustering should find 4 or 5 clusters, the 3 highest frequency clusters are
freqs = [0.362204, 0.202305, 0.111506]

srand(1)

data = readcsv("data.csv", header = true)
y = data[1][:, 1]
N = data[1][:, 2]

@time out = dpclustgibbs(y, N, iterations = 20000)

@test isapprox(sort(out.clonefrequencies, rev = true)[1:3], [0.362204, 0.202305, 0.111506], atol = 0.01)

#check no errors from plotting
plotresults(out; save = true)
@test isfile("DPclustering.png")
rm("DPclustering.png")

#check show function works
show(out)
