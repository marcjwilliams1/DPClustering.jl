# DPClustering

[![Build Status](https://travis-ci.org/marcjwilliams1/DPClustering.jl.svg?branch=master)](https://travis-ci.org/marcjwilliams1/DPClustering.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/marcjwilliams1/DPClustering.jl?branch=master&svg=true)](https://ci.appveyor.com/project/marcjwilliams1/DPClustering-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/marcjwilliams1/DPClustering.jl/badge.svg?branch=master)](https://coveralls.io/github/marcjwilliams1/DPClustering.jl?branch=master)
[![codecov.io](http://codecov.io/github/marcjwilliams1/DPClustering.jl/coverage.svg?branch=master)](http://codecov.io/github/marcjwilliams1/DPClustering.jl?branch=master)

Perform Dirichlet clustering on Varaint Allele Frequncies (VAFs) from sequencing data of cancers a la Nik-Zainal et al.

## Getting started
Package is written in the [Julia](https://julialang.org/) programming language.

To download this package use the ```Pkg.add``` function as below, which will download the package and install all the dependencies.
```
Pkg.add("DPClustering")
```

## Clustering
Clustering is invoked using the ```dpclustering``` function which takes 2 vectors of equal size: y - the number of reads reporting the mutation and N - the depth of coverage at that locus. With these, clustering can be performed and the function will return a ```DPresults``` type. There are a number of optional arguments which are all set to reasonable defaults, you may want to change the number of iterations or set verbose to false. the default is to show a log of the time taken in the gibbs sampler. To see the optional arguments and their defaults use ```?dpclustgibbs``` in the julia repl.

```
dp = dpclustering(y, N, iterations = 10000, verbose = false);
```

You can then summarise and plot the output using ```show(dp)``` and ```plotresults(dp)```.

At the moment, clustering will only work with single samples and mutations in copy neutral regions. So the input mutations should either be filtered for copy number alterations or corrected for copy number before inputting.

## Example
There is some example data provided originally in Nik-Zainal et al in the examples folder. So an analysis would proceed as follows. We'll use the Gadfly package to save a plot.
```
using DPClustering
data = readcsv("example/data.csv", header = true)
y = data[1][:, 1]
N = data[1][:, 2]

out = dpclustering(y, N)
show(out)
myplot = plotresults(out, save = true)
```

![plot](/example/example.png)

## Speed
Due to Julia's just in time compilation the sampling is relatively fast. For example, the analysis above (600 mutations) the time taken to generate 10,000 iterations/samples should be on the order of 2-3 minutes on a reasonably specced laptop.

## Acknowledgments
The model used in the Gibbs sampler is as described in Nik-Zainal et al. Bugs code provided in the supplementary information of this publication was taken as inspiration, as was code available from David Wedge's group (https://github.com/Wedge-Oxford/dpclust_docker).
