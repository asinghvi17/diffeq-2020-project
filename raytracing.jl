scene, layout = layoutscene(5)

ax = layout[1, 1] = LAxis(scene; xlabel = "𝑧", ylabel = "𝑟")

function r(z; a = 1, b = 0, c = 1)
    return if z < 0
        a
    else
        a * cos(c * z) + b / c * sin(c * z)
    end
end
zs = LinRange(-2, 10, 1000)
lines!(zs, r.(zs))

yvals = ((x, y) -> y).(LinRange(-.99, .99, 1000), LinRange(-.99, .99, 1000)')
c = 1
zs = @. 1 / (sqrt(1 - c * abs(yvals)))
heatmap!(-1..1, -1..1, zs; interpolate = true, colormap = cgrad(:viridis; alpha = .8))
