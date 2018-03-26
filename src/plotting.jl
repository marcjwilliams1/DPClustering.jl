"""
    plotresults(dp; <keyword arguments>)

Plot results from DPclustering object. Will plot histogram of raw data with density estimates from Gibbs sampling.
...
## Arguments
- `save = false`: Set to `true` if you want the plot to be saved
- `dir = ""`: Directory where the plot will be saved to. Default is the current working directory.
- `plotname = "DPclustering"`: Name to call plot when saving.
...
"""
function plotresults(dp; save = false, dir = "", plotname = "DPclustering")

  DFres = dp.DF
  DF = dp.data.DF

  l1 = layer(DFres, x = :x, y = :mean, ymin = :lq, ymax = :uq, Geom.line, Geom.ribbon,
  Theme(default_color = RGBA(0.325, 0.525, 0.608),
  lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.5)))
  l2 = layer(DF, x = :VAF, Geom.histogram(bincount=100, density = true),
  Theme(default_color = RGBA(0.5, 0.5, 0.5, 0.8)))

  myplot = plot(l1, l2,
  Guide.xlabel("VAF"),
  Guide.ylabel("Density"))

  if save == true
      Gadfly.draw(PNG(joinpath(dir, "$(plotname).png"), 4inch, 3inch), myplot)
  end

  return myplot
end
