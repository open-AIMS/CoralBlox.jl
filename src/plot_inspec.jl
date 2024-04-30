using Plots, Colors

"""
This function will plot one figure for each timestep
"""
function plot_size_class(size_classes::Matrix{SizeClass}, timestep::Int64, coral_spec)
    plots = Array{Plots.Plot}(undef, 6)
    for func_group in 1:6
        colors = sequential_palette(func_group * 360 / 6, 7)[2:end]
        plots[func_group] = plot(title=coral_spec.f_groups[func_group], ylims=(1, 1e7))
        n_corals::Float64 = 0.0
        for (size, sc) in enumerate(size_classes[func_group, :])
            _x = vcat([[cv.interval[1], cv.interval[2]] for cv in sc.cover_blocks]...)
            _y = vcat([fill(cv.diameter_density * (Δinterval(cv.interval)), 2) for cv in sc.cover_blocks]...)
            n_corals += sum(_y) / 2
            plot!(
                _x,
                _y,
                yscale=:log10,
                color=colors[size],
                fill=(0, 1, colors[size]),
                label="Size: $size",
                leg_title="SizeClass"
            )
        end

        plot!(
            [coral_spec.bins[1], coral_spec.bins[end]],
            [n_corals, n_corals],
            yscale=:log10,
            line=(2, :solid, :black)
        )
        xlabel!("Diameter (m)")
        ylabel!("Number of Corals per Size Class")
    end

    p = plot(plots..., layout=(3, 2), size=(1200, 1200), plot_title="timestep = $timestep")
    savefig(p, "./figures/gif/t$(timestep).png")
    return nothing
end
