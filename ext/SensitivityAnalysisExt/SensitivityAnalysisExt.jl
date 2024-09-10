module SensitivityAnalysisExt
using DataFrames
using Distributed
using Distributions
using ProgressMeter
using Zarr
using CoralBlox
using CoralBlox: FunctionalGroup

"""
# Functional groups simulated
- tabular_Acropora,
- corymbose_Acropora,
- corymbose_non_Acropora,  # and Pocillopora
- small_massives,
- large_massives

# Params
- scale_threshold
Threshold above which we apply the linear extension scale factor

- size_class_weight_fg_sc           (1 < fg < 5; 1 < sc < 6)
Weight for each size class size. The maximum size is fixed. This parameter only controls the
size of each size class in relation to the others (from the same functional group).
The terminal size class is ignored since it has a fixed size (0.1).

- initial_cover_weight_fg           (1 < fg < 5)
Weight for each functional group initial cover. The absolute initial cover is fixed. This
param controls the proportion of that cover that goes to each functional group.

- linear_extension_perturb_fg_sc    (1 < fg < 5; 1 < sc < 7)
Linear extension perturbation (as a fraction of the base value).

- survival_rate_perturb_fg_sc       (1 < fg < 5; 1 < sc < 7)
Survival rate perturbation.

- fecundity_weights_fg              (1 < fg < 5)
Weight params used to simulate a dummy coral reproduction scheme.
"""

# In meters
_functional_group_max_sizes()::Vector{Float64} = [1.5, 1.0, 0.5, 1.0, 1.0]

# TODO Incorporate n_size_classes in the SA
#? How?
_n_size_classes()::Int64 = 7

"""
Base linear extensions in m
"""
_base_linear_extensions() = [
    0.609456 1.07184 2.55149 5.07988 9.45091 16.8505 0.0;           # Tabular Acropora
    0.768556 1.22085 1.86447 2.82297 3.52938 3.00422 0.0;           # Corymbose Acropora
    0.190455 0.343747 0.615467 0.97477 1.70079 2.91729 0.0;         # Corymbose non-Acropora
    0.318034 0.47385 0.683729 0.710587 0.581085 0.581085 0.0;       # Small massives
    0.122478 0.217702 0.382098 0.718781 1.24172 2.08546 0.0         # Large massives
] .* 0.01

# TODO explain this -0.1 that is related to the perturbation
_base_survival_rate() = [
    0.6 0.76 0.805 0.76 0.85 0.86 0.86;                             # Tabular Acropora
    0.6 0.76 0.77 0.875 0.83 0.90 0.90;                             # Corymbose Acropora
    0.52 0.77 0.77 0.875 0.89 0.97621179 0.97621179;                # Corymbose non-Acropora
    0.72 0.87 0.77 0.98 0.996931548 0.996931548 0.996931548;        # Small massives
    0.58 0.87 0.78 0.983568572 0.984667677 0.984667677 0.984667677  # Large massives
] .- 0.1

"""
    _size_class_bounds(size_class_weights::DataFrame)

Size class bounds. Rows: functional groups; Cols: size classes.
"""
function _size_class_bounds(size_class_weights::DataFrame)::Matrix{Float64}
    max_sizes::Vector{Float64} = _functional_group_max_sizes()
    # n_size_classes == length(size_class_weights) + 1 because the last size class
    # is added by hand and it has fixed size == 0.1 and linear extension == 0.0
    n_size_classes = _n_size_classes()

    # The number of bounds is n_size_classes + 1
    size_class_bounds::Matrix{Float64} = zeros(length(max_sizes), n_size_classes + 1)
    for (idx, max_size) in enumerate(max_sizes)
        # Select all size_class_weights for current functional group
        start_idx = (1 + (idx - 1) * (n_size_classes - 1))
        end_idx = start_idx + n_size_classes - 2
        weights::Vector{Float64} = collect(size_class_weights[!, start_idx:end_idx][1, :])

        # Normalize weights using max_size
        sum_weights::Float64 = sum(weights)
        normalized_weights::Vector{Float64} = [(weight / sum_weights) * max_size for weight in weights]

        # Fill size_class_bounds. First bound is always 0.0
        size_class_bounds[idx, 2:end-1] = round.(cumsum(normalized_weights), digits=2)
        size_class_bounds[idx, end] = size_class_bounds[idx, end-1] + 0.1
    end

    return size_class_bounds
end

function _initial_cover!(
    C_cover_1::SubArray{Float64,2},
    size_class_bounds::Matrix{Float64},
    asbolute_initial_cover::Float64,
    initial_cover_weights::DataFrame
)::Matrix{Float64}
    # Normalize initial cover weights
    icc_weights::Vector{Float64} = Vector(initial_cover_weights[1, :])
    icc_weights ./= sum(icc_weights)
    functional_groups_initial_cover = asbolute_initial_cover .* icc_weights

    # Cover Weights
    n_fg::Int64 = size(size_class_bounds, 1)
    size_dists::Vector{Exponential{Float64}} = Exponential.(fill(0.5, n_fg))
    cover_weights = cdf.(size_dists, size_class_bounds[:, 2:end]) .- cdf.(size_dists, size_class_bounds[:, 1:end-1])

    #  Normalize cover_weights
    cover_weights = cover_weights ./ sum(cover_weights, dims=2)

    # Set initial cover
    C_cover_1 .= hcat(((functional_groups_initial_cover) .* eachrow(cover_weights))...)'
end

function _linear_extensions(linear_extension_perturb::DataFrame)::Matrix{Float64}
    base_linear_extensions = _base_linear_extensions()

    n_functional_groups::Int64, n_size_classes::Int64 = size(base_linear_extensions)

    return base_linear_extensions .+ (
        base_linear_extensions .*
        reshape(
            collect(linear_extension_perturb[1, :]),
            n_functional_groups,
            n_size_classes
        )
    )
end

# TODO Think about how to do this
function _survival_rates(survival_rate_perturb::DataFrame)::Matrix{Float64}
    n_functional_groups::Int64, n_size_classes::Int64 = size(_base_survival_rate())

    return _base_survival_rate() .+ (
        reshape(
        collect(survival_rate_perturb[1, :]),
        n_functional_groups,
        n_size_classes
    )
    )
end

function _run_scenario(
    params::DataFrame,
    C_cover::SubArray{Float64,3},
    habitable_area::Float64,
    asbolute_initial_cover::Float64
)
    n_timesteps::Int64, n_functional_groups::Int64, n_size_classes::Int64 = size(C_cover)

    # Each row contains SizeClass' bounds in m
    # TODO Cache this matrix
    size_class_bounds::Matrix{Float64} = _size_class_bounds(params[!, r"size_class_weight"])

    # Setup initial coral cover
    initial_cover_weights::DataFrame = params[!, r"initial_cover_weight"]
    # Fill first timestep with initial cover
    _initial_cover!(
        @view(C_cover[1, :, :]), size_class_bounds, asbolute_initial_cover, initial_cover_weights
    )
    # C_cover[2:end, :, :] .= 0.0

    # Setup linear extensions
    linear_extensions::Matrix{Float64} = _linear_extensions(
        params[!, r"linear_extension_perturb"]
    )

    # Seup survival rate
    survival_rates::Matrix{Float64} = _survival_rates(params[!, r"survival_rate_perturb"])

    habitable_max_projected_cover::Float64 = CoralBlox.max_projected_cover(
        linear_extensions,
        size_class_bounds,
        habitable_area
    )

    fecundity_weights::Vector{Float64} = collect(params[!, r"fecundity_weights"][1, :])

    # Create functional groups
    functional_groups::Vector{FunctionalGroup} = CoralBlox.FunctionalGroup.(
        eachrow(size_class_bounds[:, 1:end-1]),      # lower bounds
        eachrow(size_class_bounds[:, 2:end]),        # upper bounds
        eachrow(C_cover[1, :, :])                # initial coral covers
    )

    scale_threshold = params.scale_threshold[1] * habitable_area
    # Linear extension scale factor
    local lin_ext_scale_factors::Float64

    for tstep::Int64 in 2:n_timesteps
        # Only re-scale if total cover is above a given threshold
        lin_ext_scale_factor = if sum(C_cover[tstep-1, :, :]) < scale_threshold
            lin_ext_scale_factor = 1
        else
            CoralBlox.linear_extension_scale_factors(
                C_cover[tstep-1, :, :],
                habitable_area,
                linear_extensions,
                size_class_bounds,
                habitable_max_projected_cover,
            )

        end

        # Use scale factor to calculate growth rate
        growth_rate::Matrix{Float64} = linear_extensions .* lin_ext_scale_factor

        # Mock recruits proportional to each functional group's cover and available space
        available_frac::Float64 = (habitable_area - sum(C_cover[tstep-1, :, :])) / habitable_area
        adults_cover::Vector{Float64} = dropdims(sum(C_cover[tstep-1, :, 2:end], dims=2), dims=2)

        recruits::Vector{Float64} = adults_cover .*
                                    fecundity_weights .*
                                    log(2, 1 + available_frac)

        # Perform timestep
        CoralBlox.timestep!(
            functional_groups,
            recruits,
            growth_rate,
            survival_rates
        )

        # Write to the cover matrix
        C_cover_t = @view(C_cover[tstep, :, :])
        @assert sum(C_cover_t) <= habitable_area
        CoralBlox.coral_cover(functional_groups, C_cover_t)
    end
end

function _save_cover(C_cover::Array{Float64,4}, result_path::String)::Nothing
    cover_size = size(C_cover)
    z1 = zcreate(
        Float64,
        cover_size...,
        path=result_path,
        chunks=(1000, cover_size[2:end]...)
    )
    z1 .= C_cover
    return nothing
end

function CoralBlox.SensitivityAnalysis.run_scenarios(
    params::DataFrame;
    n_timesteps::Int64=2,
    result_path::String="cbsa_cover.zarr"
)::Array{Float64,4}
    n_scenarios::Int64 = nrow(params)
    n_functional_groups::Int64 = length(_functional_group_max_sizes())
    n_size_classes::Int64 = _n_size_classes()

    # Coral cover cache
    C_cover::Array{Float64,4} = zeros(n_scenarios, n_timesteps, n_functional_groups, n_size_classes)

    # Habitable area is the maximum possible cover that a location can hold
    habitable_area::Float64 = 1e6

    # Absolute Coral Cover is a fixed fraction of the habitable_area.
    # In the future I can vary this as part of the SA as well
    asbolute_initial_cover::Float64 = habitable_area * 0.5

    @info "Running $n_scenarios scenarios..."

    p_run_scenario(x) = _run_scenario(
        DataFrame(params[x, :]),
        @view(C_cover[x, :, :, :]),
        habitable_area,
        asbolute_initial_cover,
    )

    @showprogress pmap(p_run_scenario, 1:size(params, 1))

    @info "Saving cover to zarr file..."
    _save_cover(C_cover, result_path)

    @info "Success!"

    return C_cover
end
end
