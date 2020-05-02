using AbstractPlotting, CairoMakie, MakieLayout
using Roots
AbstractPlotting.inline!(true)

AbstractPlotting.set_theme!(
    resolution = (155, 155),
    textsize = 10,
    fontsize = 10,
    linewidth = 1,
    font = "CMU Serif Roman",
    LAxis = (
        spinewidth = 0.5,
        xgridwidth = 0.4,
        ygridwidth = 0.4,
        xticksize = -3,
        yticksize = -3,
        titlegap  = 5.
    )
)


function cycloid(xf, yf; N = 1000, g = 9.81)

    # First, we find the final angle θ₂ by
    # optimizing a loss function numerically.
    f(θ) = yf/xf - (1 - cos(θ))/(θ - sin(θ))
    θ₂ = Roots.fzero(f, π/2)
    @show θ₂

    R = yf / (1 - cos(θ₂))

    θs = LinRange(0, θ₂, N)

    return Point2f0.(
        R .* (θs .- sin.(θs)),
        R .* (1  .- cos.(θs))
    )
end

scene, layout = layoutscene(10)

ax = layout[1, 1] = LAxis(scene)
ax.title = "Brachistochrone trajectory"
ax.yreversed = true

lines!(ax, cycloid(1, .65))
# plot the straight-line path and hte
lines!(ax, [Point2f0(0), Point2f0(0, .65), Point2f0(1, .65)]; linestyle = :dash, linewidth = .7)
lines!(ax, [Point2f0(0), Point2f0(1, .65)]; linestyle = :dash, linewidth = .7)
save("brachistochrone.pdf", scene)
