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


function r(z; a = .6, b = 0, c = .1)
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

################################################################################
#                               GRIN index ray                                 #
################################################################################

function ray(u, p, t)
    r,ṙ = u

    return @SVector[ṙ, -p^2 * r / (1-p^2*r^2)]*(1-ṙ)
end

function r_analytical(z; a = 0.6, b = 0, c = 0.1)
    return if z < 0
        b*z+a
    else
        a * cos(c * z) + b / c * sin(c * z)
    end
end

scene, layout = layoutscene(3, resolution = (340, 180))
ax = layout[1, 1] = LAxis(
    scene;
    title = "Straight light rays",
    # ax.title = "Angled Light"
    xlabel = "𝑧",
    ylabel = "𝑟(𝑧)",
    xlabelpadding = .5,
    ylabelpadding = .5,
    xticklabelpad = 2,
    yticklabelpad = 2,
    xticklabelsize = 12,
    yticklabelsize = 12,
)
ax2 = layout[1, 2] = LAxis(
    scene;
    width = Fixed(40),
)
hidedecorations!(ax2)
ax2.xgridvisible = true
ax2.xticklabelsvisible = false
ax2.ygridvisible = true
ax2.yticklabelsvisible = false
ax2.title = "𝜂(𝑟)"
ax2.yaxisposition = :right
ax2.xlabel = "𝑟"
ax2.ylabel = "𝜂"
ri_range = LinRange(-0.6, 0.6, 1000)
lines!(ax2, (y -> 1/sqrt(1-(.1)^2 * y^2)).(ri_range), ri_range)


save("a.png", scene; px_per_unit = 10)

zs = LinRange(-5, 100, 1000)
lines!(ax, zs, r.(zs); linewidth = 1, color = AbstractPlotting.wong_colors[2])
lines!(ax, zs, r.(zs; a = -0.1); linewidth = 1, color = AbstractPlotting.wong_colors[1])
lines!(ax, zs, r.(zs; a = -0.4); linewidth = 1, color = AbstractPlotting.wong_colors[3])
# lines!(ax, zs, r.(zs; a = -0.6, b = 0.1); linewidth = 1, color = AbstractPlotting.wong_colors[2])
# lines!(ax, zs, r.(zs; a = 0.4, b = -0.2); linewidth = 1, color = AbstractPlotting.wong_colors[2])


function rn(z; a = 0.6, b = 0.0, c = 0.1)
    prob = ODEProblem(ray, [a, b], (0.0, 100.0), p = c)
    sol = solve(prob; save_idxs = 1)
    return (i -> if i < 0
        return b*i+a
    else
        return sol(i)
    end).(z)
end
# lines!(ax, zs, rn.(zs); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
# lines!(ax, zs, rn.(zs; a = -0.1); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
# lines!(ax, zs, rn.(zs; a = -0.4); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
# lines!(ax, zs, rn.(zs; a = -0.6, b = 0.1); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
# lines!(ax, zs, rn.(zs; a = 0.4, b = -0.2); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])

save("ray_comp.pdf", scene)
# save("ray_comp_angled.pdf", scene)


scene_angled, layout_angled = layoutscene(3, resolution = (340, 180))
ax_angled = layout_angled[1, 1] = LAxis(
    scene_angled;
    title = "Angled light rays",
    xlabel = "𝑧",
    ylabel = "𝑟(𝑧)",
    xlabelpadding = .5,
    ylabelpadding = .5,
    xticklabelpad = 2,
    yticklabelpad = 2,
    xticklabelsize = 12,
    yticklabelsize = 12,
)
ax2_angled = layout_angled[1, 2] = LAxis(
    scene_angled;
    width = Fixed(40),
)
hidedecorations!(ax2_angled)
ax2_angled.xgridvisible = true
ax2_angled.xticklabelsvisible = false
ax2_angled.ygridvisible = true
ax2_angled.yticklabelsvisible = false
ax2_angled.title = "𝜂(𝑟)"
ax2_angled.yaxisposition = :right
ax2_angled.xlabel = "𝑟"
ax2_angled.ylabel = "𝜂"
ri_angled_range = LinRange(-2, 2, 1000)
lines!(ax2_angled, (y -> 1/sqrt(1-(.1)^2 * y^2)).(ri_angled_range), ri_angled_range)

lines!(ax_angled, zs, rn.(zs; a = -0.6, b = 0.1); linewidth = 1, color = AbstractPlotting.wong_colors[2])
lines!(ax_angled, zs, rn.(zs; a = 0.4, b = -0.2); linewidth = 1, color = AbstractPlotting.wong_colors[1])

save("ray_comp_angled.pdf", scene_angled)

save("ray_comp_angled.png", scene_angled; px_per_unit = 10)


## anal sim comp

scene_comp, layout_comp = layoutscene(3, resolution = (400, 250))
ax_comp = layout_comp[1, 1] = LAxis(scene_comp)
ax_comp.title = "Straight Light"
# ax_comp.title = "Angled Light"
ax_comp.xlabel = "𝑧"
ax_comp.ylabel = "𝑟(𝑧)"
ax_comp.xlabelpadding = .5
ax_comp.ylabelpadding = .5
ax_comp.xticklabelpad = 2
ax_comp.yticklabelpad = 2
ax_comp.xticklabelsize = 12
ax_comp.yticklabelsize = 12

function r(z; a = 0.6, b = 0, c = 0.1)
    return if z < 0
        b*z+a
    else
        a * cos(c * z) + b / c * sin(c * z)
    end
end
zs = LinRange(-5, 100, 1000)
lines!(ax_comp, zs, r.(zs); linewidth = 1, color = AbstractPlotting.wong_colors[2])
lines!(ax_comp, zs, r.(zs; a = -0.1); linewidth = 1, color = AbstractPlotting.wong_colors[2])
lines!(ax_comp, zs, r.(zs; a = -0.4); linewidth = 1, color = AbstractPlotting.wong_colors[2])
# lines!(ax_comp, zs, r.(zs; a = -0.6, b = 0.1); linewidth = 1, color = AbstractPlotting.wong_colors[2])
# lines!(ax_comp, zs, r.(zs; a = 0.4, b = -0.2); linewidth = 1, color = AbstractPlotting.wong_colors[2])


function rn(z; a = 0.6, b = 0.0, c = 0.1)
    if z < 0
        return b*z+a
    else
        prob = ODEProblem(ray, [a, b], (0.0, 100.0), p = c)
        sol = solve(prob; save_idxs = 1)
        return sol(z)
    end
end
lines!(ax_comp, zs, rn.(zs); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
lines!(ax_comp, zs, rn.(zs; a = -0.1); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
lines!(ax_comp, zs, rn.(zs; a = -0.4); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
# lines!(ax_comp, zs, rn.(zs; a = -0.6, b = 0.1); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
# lines!(ax_comp, zs, rn.(zs; a = 0.4, b = -0.2); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])

save("ray_comp.pdf", scene_comp)



scene_comp_angled, layout_comp_angled = layoutscene(3, resolution = (400, 250))
ax_comp_angled = layout_comp_angled[1, 1] = LAxis(scene_comp_angled)
ax_comp_angled.title = "Straight Light"
# ax_comp_angled.title = "Angled Light"
ax_comp_angled.xlabel = "𝑧"
ax_comp_angled.ylabel = "𝑟(𝑧)"
ax_comp_angled.xlabelpadding = .5
ax_comp_angled.ylabelpadding = .5
ax_comp_angled.xticklabelpad = 2
ax_comp_angled.yticklabelpad = 2
ax_comp_angled.xticklabelsize = 12
ax_comp_angled.yticklabelsize = 12

function r(z; a = 0.6, b = 0, c = 0.1)
    return if z < 0
        b*z+a
    else
        a * cos(c * z) + b / c * sin(c * z)
    end
end
zs = LinRange(-5, 100, 1000)
# lines!(ax_comp_angled, zs, r.(zs); linewidth = 1, color = AbstractPlotting.wong_colors[2])
# lines!(ax_comp_angled, zs, r.(zs; a = -0.1); linewidth = 1, color = AbstractPlotting.wong_colors[2])
# lines!(ax_comp_angled, zs, r.(zs; a = -0.4); linewidth = 1, color = AbstractPlotting.wong_colors[2])
lines!(ax_comp_angled, zs, r.(zs; a = -0.6, b = 0.1); linewidth = 1, color = AbstractPlotting.wong_colors[2])
lines!(ax_comp_angled, zs, r.(zs; a = 0.4, b = -0.2); linewidth = 1, color = AbstractPlotting.wong_colors[2])


function rn(z; a = 0.6, b = 0.0, c = 0.1)
    if z < 0
        return b*z+a
    else
        prob = ODEProblem(ray, [a, b], (0.0, 100.0), p = c)
        sol = solve(prob; save_idxs = 1)
        return sol(z)
    end
end
# lines!(ax_comp_angled, zs, rn.(zs); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
# lines!(ax_comp_angled, zs, rn.(zs; a = -0.1); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
# lines!(ax_comp_angled, zs, rn.(zs; a = -0.4); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
lines!(ax_comp_angled, zs, rn.(zs; a = -0.6, b = 0.1); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])
lines!(ax_comp_angled, zs, rn.(zs; a = 0.4, b = -0.2); linewidth = 0.6, color = AbstractPlotting.wong_colors[1])

save("ray_comp.pdf", scene_comp)
