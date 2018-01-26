#Running DP clustering should find 4 clusters at the following frequencies:
freqs = [0.0389717, 0.0475582, 0.111506, 0.202305, 0.362204]

srand(1)

data = readcsv("data.csv", header = true)
y = data[1][:, 1]
N = data[1][:, 2]

out = dpclustgibbs(y, N, iterations = 20000)

@test isapprox(out.clonefrequencies, freqs, atol = 0.01)
