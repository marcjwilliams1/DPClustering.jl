function plotresults(dp)

  DFres = dp.DF
  DF = dp.data.DF

  l1 = layer(DFres, x = :x, y = :mean, ymin = :lq, ymax = :uq, Geom.line, Geom.ribbon,
  Theme(default_color = RGBA(0.325, 0.525, 0.608),
  lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.5)))
  l2 = layer(DF, x = :VAF, Geom.histogram(bincount=100, density = true),
  Theme(default_color = RGBA(0.5, 0.5, 0.5, 0.8)), bar_spacing = -0.05)

  myplot = plot(l1, l2,
  Guide.xlabel("VAF"),
  Guide.ylabel("Density"))

  return myplot
end
