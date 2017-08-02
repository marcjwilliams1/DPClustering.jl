# DPclustering

Perform Dirichlet clustering on Varaint Allele Frequncies (VAFs) from sequencing data of cancers a la Nik-Zainal et al.

## Getting started
Package is written in the [Julia](https://julialang.org/) programming language, and has been tested with julia v0.5.1.

To download this package use the ```Pkg.clone``` function as below, which will download the package and install all the dependencies.
```
Pkg.clone("https://github.com/marcjwilliams1/DPclustering.jl")
```
To correctly calculate the probability density requires a weighted kernel density estimator which is only available in the latest version of ```KernelDensity.jl```. To update the version to the latest one use the following command once you have cloned the ```DPclustering``` package.
```
Pkg.checkout("KernelDensity")
```

## Clustering
Clustering is invoked using the ```dpclustgibbs``` function which takes 2 vectors of equal size: y - the number of reads reporting the mutation and N - the depth of coverage at that locus. With these, clustering can be performed and the function will return a ```DPresults``` type. There are a number of optional arguments which are all set to reasonable defaults, you may want to change the number of iterations or set verbose to false. the default is to show a log of the time taken in the gibbs sampler. To see the optional arguments and their defaults use ```?dpclustgibbs``` in the julia repl.

```
dp = dpclustgibbs(y, N, iterations = 2000, verbose = false);
```

You can then summarise and plot the output using ```show(dp)``` and ```plotresults(dp)```.

At the moment, clustering will only work with single samples and mutations in copy neutral regions. So the input mutations should either be filtered for copy number alterations or corrected for copy number before inputting.

## Example
There is some example data provided originally in Nik-Zainal et al in the examples folder. So an analysis would proceed as follows. We'll use the Gadfly package to save a plot.
```
using DPclustering
using Gadfly
data = readcsv("example/data.csv", header = true)
y = data[1][:, 1]
N = data[1][:, 2]

out = dpclustgibbs(y, N)
show(out)
myplot = plotresults(out)
draw(PNG("example/example.png", 15cm, 10cm), myplot)
```

![plot](/example/example.png)


## Acknowledgments
The model used in the Gibbs sampler is as described in Nik-Zainal et al. Bugs code provided in the supplementary information of this publication was taken as inspiration, as was code available from the Peter van Loo and David Wedge groups (eg: https://github.com/Wedge-Oxford/dpclust_docker).
