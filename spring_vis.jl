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

################################################################################
#                           3-spring, 2-mass system                            #
################################################################################

@inbounds function springmass(u, p, t)
    k₁, k₂, k₃, m₁, m₂ = p
    x₁, ẋ₁, x₂, ẋ₂ = u

    return @SVector[
        # D(x₁) = ẋ₁
        ẋ₁,
        # D(D(x₁)) = (k₁x₁ - k₂(x₂ - x₁))/m₁
        -(k₁ * x₁ + k₂ * (x₁ - x₂))/m₁,
        # D(x₂) = ẋ₂
        ẋ₂,
        # D(D(x₂)) = (k₃x₂ - k₂(x₂ - x₁))/m₂
        -(k₃ * x₂ + k₂ * (x₂ - x₁))/m₂,
    ]
end

function energy(resid, u, p, t)
    k₁, k₂, k₃, m₁, m₂ = p
    x₁, ẋ₁, x₂, ẋ₂ = u
    # This is the Lagrangian, which must be kept at 0.
    resid[1] = (m₁ * ẋ₁^2 + m₂ * ẋ₂^2)/2 - (k₁ * x₁^2 + k₃ * x₂^2 - k₂ * (x₂ - x₁)^2)
    # These are irrelevant - we only want to conserve the one Lagrangian.
    resid[2] = 0
    resid[3] = 0
    resid[4] = 0
end


trange = (0.0, 80.0)
prob = ODEProblem(
    springmass,
    [0.0, -1.0, 0.0, -5.0],
    trange,
    p = [40.0, 10.0, 40.0, 1.0, 1.0],
    # callback = ManifoldProjection(energy)
)

# Parameter sets:
# The integrator bugs out at
# p = [40.0, 10.0, 40.0, 1.0, 1.0],
# at least at the initial stages.
sol = solve(prob, Tsit5(); save_idxs = [1, 2])

# import Plots
#
Plots.plot(sol; vars = (1, 2))
Plots.plot(sol; vars = [1, 2])




scene, layout = layoutscene(5, resolution = (200, 170), font = "CMU Serif Roman")
ax = layout[1, 1]  = LAxis(scene)
ax.xlabel = "𝑡"
ax.ylabel = "𝑥"
ax.xticklabelpad = 2
ax.yticklabelpad = 2
ax.xlabelpadding = 0
ax.ylabelpadding = 0

lines!(ax, sol.t, sol[1, :]; color = AbstractPlotting.wong_colors[2], linewidth = .5)
lines!(ax, sol.t, sol[2, :]; color = AbstractPlotting.wong_colors[1], linewidth = .5)

leg = layout[2, 1] = LLegend(
    scene,
    ax.scene.plots,
    ["𝑥1", "𝑥2"];
    tellwidth = false, rowgap = 0.0,
    framewidth = 0, padding = (0,0,0,0),
    patchlabelgap = 3, patchsize = (5, 5),
    linewidth = 1
)

rowsize!(layout, 2, Fixed(20))
rowgap!(layout, 0)
save("3spring2mass.pdf", scene)
