using CairoMakie, AbstractPlotting, MakieLayout
using DifferentialEquations, StaticArrays

function spring(u, p, t)
    m, k = p
    return @SVector[
        u[2],
        -k/m * u[1]
    ]
end

prob = ODEProblem(spring, [-1.0, 0.0], (0.0, 2.0), p = [1.0, 490.0])

sol = solve(prob; save_idxs = 1)

ts = LinRange(extrema(sol.t)..., 3000)
xs = sol.(ts)
scene, layout = layoutscene(5)
ax = layout[1, 1] = LAxis(scene)
ax.title = "Spring-mass timeseries"
ax.xlabel = "Time"
ax.ylabel = "𝑥(𝑡)"
ax.xlabelpadding = .5
ax.ylabelpadding = .5
ax.xticklabelpad = 1
ax.yticklabelpad = 1
ax.xticklabelsize = 8
ax.yticklabelsize = 8
lines!(ax, ts, xs; linewidth = .5)

save("spring_mass_1.pdf", scene)
