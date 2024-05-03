using Plots, Colors

"""
This function will plot one figure for each timestep
"""
function plot_size_class(size_classes::Matrix{SizeClass}, timestep::Int64, coral_spec)
    n_species = length(coral_spec.f_groups)
    n_bins = size(coral_spec.bins)[2] - 1
    plots = Array{Plots.Plot}(undef, n_species)
    for species in 1:n_species
        colors = sequential_palette(species * 360 / n_species, n_bins + 1)[2:end]
        plots[species] = plot(title=coral_spec.f_groups[species], ylims=(1, 1e7))
        n_corals::Float64 = 0.0

        for (size, sc) in enumerate(size_classes[species, :])
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
