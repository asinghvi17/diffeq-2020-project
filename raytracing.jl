scene, layout = layoutscene(3; resolution = (200, 100))

ax = layout[1, 1] = LAxis(scene)
ax.xlabel = "𝑧"
ax.ylabel = "𝑟"
ax.xlabelpadding = .5
ax.ylabelpadding = .5
ax.xticklabelpad = 2
ax.yticklabelpad = 2
ax.xticklabelsize = 8
ax.yticklabelsize = 8

# show the cylinder's refractive index
yvals = ((x, y) -> y).(LinRange(-.99, .99, 100), LinRange(-.99, .99, 100)')
c = 1
colors = @. 1 / (sqrt(1 - c * abs(yvals)))
heatmap!(ax, 0..10, -1..1, colors; interpolate = true, colormap = cgrad(:grays; alpha = .5))


function r(z; a = 1, b = 0, c = 1)
    return if z < 0
        a
    else
        a * cos(c * z) + b / c * sin(c * z)
    end
end
zs = LinRange(-2, 10, 1000)
lines!(ax, zs, r.(zs); linewidth = .5, color = AbstractPlotting.wong_colors[2])
lines!(ax, zs, r.(zs; a = -.3); linewidth = .5, color = AbstractPlotting.wong_colors[1])

save("rays.png", scene; px_per_unit = 10)
save("rays.eps", scene)
save("rays.svg", scene)
save("rays.pdf", scene)
