# The runnable Quick Start example referenced from README.md. Generates docs/assets/plain_cb.png.
# Run standalone with `julia --project=docs docs/make_quickstart_plot.jl` from the repo root.

using PlotlyLight
using PlotlyKaleido
using EasyConfig: Config
using CoralBlox: FunctionalGroup
using CoralBlox: linear_extension_scale_factors, max_projected_cover
using CoralBlox: timestep!, coral_cover

"""
Runs the Quick Start example for a given `scale_threshold` (as a fraction of
`habitable_area`) and returns the relative cover per functional group and the total
relative cover over time.
"""
function run_quickstart(scale_threshold_frac::Float64)
    n_timesteps::Int64 = 100
    n_functional_groups::Int64 = 3
    n_size_classes::Int64 = 4

    # Coral cover cache
    C_cover::Array{Float64,3} = zeros(n_timesteps, n_functional_groups, n_size_classes)

    # Habitable area is the maximum possible cover
    habitable_area::Float64 = 1e6

    # Each row contains SizeClass' bounds
    size_class_bounds::Matrix{Float64} = [
        0 0.05 0.8 1.4 1.5
        0 0.05 0.5 0.9 1.0
        0 0.05 0.5 0.9 1.0
    ]

    # Mock initial coral cover
    C_cover[1, :, :] = [
        0.08 0.05 0.02 0.005
        0.07 0.06 0.03 0.007
        0.05 0.05 0.02 0.003
    ] .* habitable_area

    # Create functional groups
    functional_groups::Vector{FunctionalGroup} = FunctionalGroup.(
        eachrow(size_class_bounds[:, 1:(end - 1)]),      # lower bounds
        eachrow(size_class_bounds[:, 2:end]),        # upper bounds
        eachrow(C_cover[1, :, :])                # initial coral covers
    )

    # Mock linear extensions
    linear_extensions::Matrix{Float64} = [
        0.004 0.025 0.1 0.0
        0.003 0.007 0.005 0.0
        0.002 0.004 0.01 0.0
    ]

    # Mock survival rate
    survival_rate::Matrix{Float64} = [
        0.5 0.6 0.65 0.8
        0.6 0.7 0.8 0.8
        0.5 0.8 0.8 0.8
    ]

    # Calculate maximum projected cover
    habitable_max_projected_cover = max_projected_cover(
        linear_extensions,
        size_class_bounds,
        habitable_area
    )

    # Only apply the spatial competition factor (linear extension scale factor) when cover is
    # above the scale_threshold. Below that threshold we assume there is enough available space
    # that spatial competition doesn't affect the dynamics, so growth is left unscaled. Note that
    # this thresholding behaviour is a choice made by this example, not a requirement of CoralBlox.
    scale_threshold = scale_threshold_frac * habitable_area

    # Spatial competition factor (linear extension scale factor)
    local scale_factor::Float64

    for tstep::Int64 in 2:n_timesteps
        # Apply the spatial competition factor to linear_extension when cover is above
        # scale_threshold, to account for corals competing for a shared, limited habitable area
        scale_factor = if sum(C_cover[tstep - 1, :, :]) < scale_threshold
            1
        else
            linear_extension_scale_factors(
                C_cover[tstep - 1, :, :],
                habitable_area,
                linear_extensions,
                size_class_bounds,
                habitable_max_projected_cover
            )
        end

        # Use the spatial competition factor to calculate growth
        growth_rate::Matrix{Float64} = linear_extensions .* scale_factor

        # Mock recruits proportional to each functional group's cover and available space
        # This mock is for example purposes only
        available_space::Float64 = habitable_area - sum(C_cover[tstep - 1, :, :])
        available_proportion::Float64 = available_space / habitable_area
        adults_cover::Vector{Float64} = dropdims(
            sum(C_cover[tstep - 1, :, 2:end], dims=2), dims=2
        )

        recruits_weights::Vector{Float64} = [0.6, 0.9, 1.5]
        availability_weight::Float64 = log(2.5, 1.5 + available_proportion)

        # In reality, recruits cover for each functional group would come from a fecundity model
        recruits::Vector{Float64} = adults_cover .* recruits_weights .* availability_weight

        # Perform timestep
        timestep!(
            functional_groups,
            recruits,
            growth_rate,
            survival_rate
        )

        # Write to the cover matrix
        coral_cover(functional_groups, @view(C_cover[tstep, :, :]))
    end

    # Relative cover per functional group over time, as a fraction of habitable area
    relative_cover = dropdims(sum(C_cover, dims=3), dims=3) ./ habitable_area
    total_relative_cover = dropdims(sum(relative_cover, dims=2), dims=2)

    return n_timesteps, n_functional_groups, relative_cover, total_relative_cover
end

# Panel comparing how the choice of scale_threshold (as a fraction of habitable_area)
# affects the space-competition growth penalty.
scale_thresholds = [0.5, 0.7, 0.9]
# Okabe-Ito colorblind-safe palette (blue, orange, bluish green); vermillion is reserved for
# the spatial competition threshold line below.
fg_colors = ["#0072B2", "#E69F00", "#009E73"]
threshold_color = "#D55E00"

p = Plot()
p.layout.title.text = "CoralBlox Quick Start example — effect of scale_threshold"
p.layout.title.y = 0.99
p.layout.title.yanchor = "top"
p.layout.title.yref = "container"
p.layout.annotations = Config[]
p.layout.margin = Config(t=90, b=45)
p.layout.width = 720
p.layout.height = 1000
p.layout.legend = Config(orientation="h", x=0.5, xanchor="center", y=1.0, yanchor="bottom")

n_panels = length(scale_thresholds)
gap = 0.03
top = 0.97  # leaves headroom above the panels for the title and horizontal legend
panel_height = (top - (n_panels - 1) * gap) / n_panels

# Single shared y-axis label, vertically centred on the panel stack
push!(p.layout.annotations, Config(
    text="Relative cover",
    xref="paper", yref="paper",
    x=-0.08, y=top / 2, showarrow=false, textangle=-90,
    font=Config(size=14)
))

for (i, thr) in enumerate(scale_thresholds)
    n_timesteps, n_functional_groups, relative_cover, total_relative_cover = run_quickstart(thr)

    suffix = i == 1 ? "" : string(i)
    xaxis_name = "x" * suffix
    yaxis_name = "y" * suffix

    # Row 1 (i == 1) goes on top, so its domain starts closest to y = top
    domain_end = top - (i - 1) * (panel_height + gap)
    domain_start = domain_end - panel_height

    for fg in 1:n_functional_groups
        p(
            x=1:n_timesteps, y=relative_cover[:, fg],
            mode="lines", name="Functional group $fg",
            legendgroup="fg$fg", showlegend=(i == 1),
            line=Config(color=fg_colors[fg]),
            xaxis=xaxis_name, yaxis=yaxis_name
        )
    end
    p(
        x=1:n_timesteps, y=total_relative_cover,
        mode="lines", name="Total cover",
        legendgroup="total", showlegend=(i == 1),
        line=Config(dash="dash", color="black"),
        xaxis=xaxis_name, yaxis=yaxis_name
    )
    p(
        x=[1, n_timesteps], y=[thr, thr],
        mode="lines", name="Spatial competition threshold",
        legendgroup="threshold", showlegend=(i == 1),
        line=Config(color=threshold_color, dash="dot"),
        xaxis=xaxis_name, yaxis=yaxis_name
    )

    xaxis = getproperty(p.layout, Symbol("xaxis" * suffix))
    xaxis.domain = [0, 1]
    xaxis.anchor = yaxis_name
    if i == n_panels
        xaxis.title.text = "Timestep"
    else
        xaxis.showticklabels = false
    end

    yaxis = getproperty(p.layout, Symbol("yaxis" * suffix))
    yaxis.domain = [domain_start, domain_end]
    yaxis.anchor = xaxis_name
    yaxis.range = [0, 1]
    yaxis.dtick = 0.1
end

PlotlyKaleido.start()
PlotlyKaleido.savefig(
    p, joinpath(@__DIR__, "assets", "plain_cb.png");
    width=p.layout.width, height=p.layout.height
)
