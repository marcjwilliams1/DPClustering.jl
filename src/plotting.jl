"""
    plotresults(dp; <keyword arguments>)

Plot results from DPClustering object. Will plot histogram of raw data with density estimates from Gibbs sampling.
...
## Arguments
- `save = false`: Set to `true` if you want the plot to be saved
- `dir = ""`: Directory where the plot will be saved to. Default is the current working directory.
- `plotname = "DPclustering"`: Name to call plot when saving.
...
"""
function plotresults(dp; save = false, dir = "", plotname = "DPClustering")

  DFres = dp.DF
  DF = dp.data.DF

  histogram(DF[:VAF], nbins = 100, normed = true, linecolor = :white,
  color = RGBA(0.42, 0.5, 0.5, 0.6))
  plot!(DFres[:x], DFres[:mean], color = RGBA(0.325, 0.525, 0.608))
  plot!(DFres[:x], DFres[:uq], fillrange = DFres[:lq],
               fillalpha = 0.5,
               fillcolor = RGBA(0.325, 0.525, 0.608), linecolor = false,
               markerstrokecolor=:white, titlefont = font(12, "Calibri"), ytickfont = font(10, "Calibri"), xtickfont = font(10, "Calibri"), legend = false, grid = false,
               yaxis = ("Density"), xaxis = ("VAF"))

  if save == true
      Plots.savefig(joinpath(dir, "$(plotname).pdf"))
  end

end
