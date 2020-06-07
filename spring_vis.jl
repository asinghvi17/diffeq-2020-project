using CairoMakie, AbstractPlotting, MakieLayout
using DifferentialEquations, StaticArrays, ModelingToolkit

function plot_spring_sol(sol)
    scene, layout = layoutscene()
    ax = layout[1, 1] = LAxis(scene)
    trange = LinRange(extrema(sol.t)..., 1000)

    data = sol(trange)

    lines!(ax, trange, data[1, :]; color = AbstractPlotting.wong_colors[2], linewidth = .5)
    lines!(ax, trange, data[2, :]; color = AbstractPlotting.wong_colors[1], linewidth = .5)

    return scene
end

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
        # D(D(x₁)) = (k₁x₁ - k₂(x₁ - x₂))/m₁
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


(t, (k₁, k₂, k₃), (m₁, m₂)) = @parameters t k[1:3] m[1:2]
(x₁, x₂), = @variables x[1:2](t)
@derivatives D'~t

eqs = [
    D(D(x₁)) ~ -(k₁ * x₁ + k₂ * (x₁ - x₂))/m₁,
    D(D(x₂)) ~ -(k₃ * x₂ + k₂ * (x₂ - x₁))/m₂,
]

sys = ode_order_lowering(ODESystem(eqs, t, [x₁, x₂], [k₁, k₂, k₃, m₁, m₂]))

u0 = [
    x₁ => -1.0,
    x₂ =>  1.0,
    D(x₁) => 5,
    D(x₂) => 0.0,
]

p = [
    m₁ => 1.0,
    m₂ => 1.0,
    k₁ => 40.0,
    k₂ => 40.0,
    k₃ => 40.0,
]

function energy(resid, u, p, t)
    m₁, m₂, k₁, k₂, k₃ = p
    x₁, ẋ₁, x₂, ẋ₂ = u
    # This is the Lagrangian, which must be kept at 0.
    resid[1] = (m₁ * ẋ₁^2 + m₂ * ẋ₂^2)/2 - (k₁ * x₁^2 + k₃ * x₂^2 - k₂ * (x₂ - x₁)^2)
    resid[2] = 0
    resid[3] = 0
    resid[4] = 0
end
tspan = (0.0, 10.0)

prob = ODEProblem(sys, u0, tspan, p; cb = ManifoldProjection(energy))

sol = solve(prob, Vern9())

plot_spring_sol(sol)


solve(prob, Vern9(); cb = ManifoldProjection(energy))
scene, layout = layoutscene(5, resolution = (200, 170), font = "CMU Serif Roman")
ax = layout[1, 1]  = LAxis(scene)
ax.xlabel = "𝑡"
ax.ylabel = "𝑥"
ax.xticklabelpad = 2
ax.yticklabelpad = 3
ax.xlabelpadding = 0
ax.ylabelpadding = 0

trange = LinRange(tspan..., 1000)
data = sol(trange)

lines!(ax, trange, data[1, :]; color = AbstractPlotting.wong_colors[2], linewidth = .5)
lines!(ax, trange, data[2, :]; color = AbstractPlotting.wong_colors[1], linewidth = .5)

leg = layout[2, 1] = LLegend(
    scene,
    ax.scene.plots,
    ["mass 1", "mass 2"];
    tellwidth = false, rowgap = 0.0,
    framewidth = 0, padding = (0,0,0,0),
    patchlabelgap = 3, patchsize = (10, 10),
    linewidth = 1,
    orientation = :horizontal,
)

rowsize!(layout, 2, Fixed(20))
rowgap!(layout, 0)
save("3spring2mass_with_velocity.pdf", scene)


## MTK experimentation
Plots.plot(sol; vars = [:x₁, :x₂])
