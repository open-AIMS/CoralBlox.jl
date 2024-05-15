module blocks_model

include("coral_spec.jl")
# TODO Use only Tuples? Maybe have an intermediate struct called CoralBlock or smth?
mutable struct CoverBlock
    diameter_density::Float64       # Number of corals / m
    interval::NTuple{2,Float64}
end

function CoverBlock(initial_cover::Float64, bin_i::Float64, bin_f::Float64)::CoverBlock
    return CoverBlock(initial_cover / area_factor(bin_i, bin_f), (bin_i, bin_f))
end

#TODO I should remove linear extension, size_class and survival_rate from here to improve performance
struct SizeClass
    cover_blocks::Vector{CoverBlock}
    interval::NTuple{2,Float64}
    size_class::Int64
    linear_extension::Float64           # m/year
    survival_rate::Float64              # %
end

"""
Builds a SizeClass. There will be only one CoverBlock (so the CoverBlock and the SizeClass
intervals will be the same)
"""
function SizeClass(
    cover_block::CoverBlock,
    size_class::Int64,
    linear_extension::Float64,
    survival_rate::Float64
)
    return SizeClass([cover_block], cover_block.interval, size_class, linear_extension, survival_rate)
end

function SizeClass(
    cover_blocks::Vector{CoverBlock},
    size_class::SizeClass
)
    return SizeClass(cover_blocks, size_class.interval, size_class.size_class, size_class.linear_extension, size_class.survival_rate)
end

""" Area factor for a corals with diameter between bin_i and bin_f """
area_factor(bin_i, bin_f) = (π / 12) * (bin_f^3 - bin_i^3)

""" Interval (or lower and upper bin bounds) of a SizeClass """
interval(size_class::SizeClass)::NTuple{2,Float64} = size_class.interval
interval(cover_block::CoverBlock)::NTuple{2,Float64} = cover_block.interval
interval_lower_bound(cover_block::CoverBlock)::Float64 = cover_block.interval[1]
interval_lower_bound(size_class::SizeClass)::Float64 = size_class.interval[1]
interval_upper_bound(cover_block::CoverBlock)::Float64 = cover_block.interval[2]
interval_upper_bound(size_class::SizeClass)::Float64 = size_class.interval[2]
Δinterval(interval::NTuple{2,Float64})::Float64 = interval[2] - interval[1]

""" Diameter Linear Extension [m/year] """
linear_extension(size_class::SizeClass)::Float64 = size_class.linear_extension

""" Settlers base fraction for each species[m²] """
#settlers_base(settlers_proportion)::Vector{Float64} = k_max() .* settlers_proportion

""" Survival rate (after background mortality) """
survival_rate(size_class::SizeClass)::Float64 = size_class.survival_rate

function size_class_cover(size_class::SizeClass)::Float64
    return sum(cover_block_cover.(size_class.cover_blocks))
end

function cover_block_cover(cover_block::CoverBlock)::Float64
    return cover_block.diameter_density * area_factor(cover_block.interval...)
end

"""
Fraction of k_max filled with each size class corals at time t
"""
function initial_cover(coral_spec)
    @assert (sum(coral_spec.initial_cover_fracs) <= 1) "Sum of all relative cover should be ≤ 1"
    return coral_spec.k_max .* coral_spec.initial_cover_fracs
end

""" Run the model :) """
function run_model(n_timesteps=100)
    coral_spec = CoralSpec()

    n_species::Int64, n_bins::Int64 = size(coral_spec.linear_extension)
    cover::Array{Float64} = zeros(Float64, n_timesteps, n_species, n_bins)
    cover[1, :, :] = initial_cover(coral_spec)

    cover_blocks::Matrix{CoverBlock} = CoverBlock.(
        cover[1, :, :],
        coral_spec.bins[:, 1:end-1],
        coral_spec.bins[:, 2:end]
    )
    size_classes::Matrix{SizeClass} = SizeClass.(
        cover_blocks,
        repeat(1:n_bins, 1, n_species)',
        coral_spec.linear_extension,
        coral_spec.survival_rate
    )

    ch = zeros(6, 7)
    apply_changes!(size_classes, @view ch[2:end, :])

    for t in 2:n_timesteps
        cover[t, :, :] .= timestep_iteration(
            t,
            cover,
            size_classes,
            coral_spec,
        )
    end

    return cover
end

""" Fraction of k_max available """
function _k_area(cover::Array{Float64}, k_max::Float64)::Float64
    sum_cov::Float64 = sum(cover)
    # sum_cov > k_max ? error("k-area should not be greater than k_max. k-area: $(sum_cov), k-max: $(k_max)") : nothing

    # Consider values close to 0.0 (e.g., 1e-214) as being 0.0
    # https://github.com/JuliaLang/julia/issues/23376#issuecomment-324649815
    if ((sum_cov - k_max) + 1.0) ≈ 1.0
        return 0.0
    elseif (sum_cov - k_max) > 0.0
        @warn "k-area should not be greater than k_max. k-area: $(sum_cov), k-max: $(k_max)"
        return 0.0
    end

    return 1.0 - (sum(cover) / k_max)
end

function new_small_size_class(
    current_size_class::SizeClass,
    recruits::Float64,
    current_growth::Float64
)::SizeClass
    new_cover_blocks = move_current_blocks(current_size_class, current_growth)
    new_cover_blocks = append!([CoverBlock(recruits, current_size_class.interval[1], current_size_class.interval[2])], new_cover_blocks)
    return SizeClass(new_cover_blocks, current_size_class)
end
function new_small_size_class(
    current_size_class::SizeClass,
    current_growth::Float64
)::SizeClass
    new_cover_blocks = move_current_blocks(current_size_class, current_growth)
    return SizeClass(new_cover_blocks, current_size_class)
end

function new_medium_size_class(
    prev_size_class::SizeClass,
    current_size_class::SizeClass,
    prev_growth::Float64,
    current_growth::Float64
)::SizeClass
    new_cover_blocks_prev = move_prev_blocks(prev_size_class, prev_growth)
    new_cover_blocks_current = move_current_blocks(current_size_class, current_growth)

    new_cover_blocks = vcat(new_cover_blocks_prev, new_cover_blocks_current)
    return SizeClass(new_cover_blocks, current_size_class)
end

function move_prev_blocks(size_class::SizeClass, growth::Float64)::Vector{CoverBlock}
    # Select blocks that are going to next SizeClass
    new_upper_bounds = interval_upper_bound.(size_class.cover_blocks) .+ growth
    sc_upper_bound = size_class.interval[2]
    blocks_within_range = new_upper_bounds .> sc_upper_bound

    new_blocks = Array{CoverBlock}(undef, sum(blocks_within_range))
    for (idx, cover_block) in enumerate(size_class.cover_blocks[blocks_within_range])
        grown_lb = interval(cover_block)[1] + growth
        new_lower_bound = grown_lb < sc_upper_bound ? sc_upper_bound : grown_lb
        new_upper_bound = interval(cover_block)[2] + growth
        new_interval = (new_lower_bound, new_upper_bound)
        new_diameter_density = cover_block.diameter_density * size_class.survival_rate
        new_blocks[idx] = CoverBlock(new_diameter_density, new_interval)
    end

    return new_blocks
end

function move_current_blocks(size_class::SizeClass, growth::Float64)::Vector{CoverBlock}
    new_lower_bounds = interval_lower_bound.(size_class.cover_blocks) .+ growth
    sc_upper_bound = size_class.interval[2]
    blocks_within_range = new_lower_bounds .< sc_upper_bound

    new_blocks = Array{CoverBlock}(undef, sum(blocks_within_range))
    for (idx, cover_block) in enumerate(size_class.cover_blocks[blocks_within_range])
        grown_ub = interval(cover_block)[2] + growth
        new_upper_bound = grown_ub < sc_upper_bound ? grown_ub : sc_upper_bound
        new_lower_bound = interval_lower_bound(cover_block) + growth
        new_interval = (new_lower_bound, new_upper_bound)
        new_diameter_density = cover_block.diameter_density * size_class.survival_rate
        new_blocks[idx] = CoverBlock(new_diameter_density, new_interval)
    end

    return new_blocks
end

function new_large_size_class(
    prev_size_class::SizeClass,
    current_size_class::SizeClass,
    prev_growth::Float64,
)
    new_upper_bound = interval_upper_bound.(prev_size_class.cover_blocks) .+ prev_growth
    blocks_within_range = new_upper_bound .> prev_size_class.interval[2]

    n_new_corals::Float64 = 0.0
    for cover_block in prev_size_class.cover_blocks[blocks_within_range]
        new_diameter_density = cover_block.diameter_density * prev_size_class.survival_rate
        interval_width = cover_block.interval[2] + prev_growth - max(
            cover_block.interval[1] + prev_growth, current_size_class.interval[1]
        )
        n_new_corals += new_diameter_density * interval_width
    end

    current_Δinterval = current_size_class.interval[2] - current_size_class.interval[1]
    current_cover_block = current_size_class.cover_blocks[1]
    new_diameter_density = (n_new_corals / current_Δinterval) + current_cover_block.diameter_density * current_size_class.survival_rate

    new_cover_block = CoverBlock(new_diameter_density, current_cover_block.interval)
    return SizeClass([new_cover_block], current_size_class)
end

function virtual_cover(size_classes, coral_spec)::Float64
    return sum(size_class_virtual_cover.(size_classes, coral_spec.linear_extension))
end

function size_class_virtual_cover(size_class, linear_extension)
    return sum(cover_block_virtual_cover.(size_class.cover_blocks, [linear_extension]))
end

function cover_block_virtual_cover(cover_block, linear_extension)
    return cover_block.diameter_density * Δinterval(cover_block.interval) * linear_extension
end

"""
- `cover` : array with dimensions timesteps X species X size class representing coral cover
- `t` : Time step of that iteration
- `bins` :
- `linear_extension` :
- `base_mortality_rate` :
"""
function timestep_iteration(
    t::Int64,
    cover::Array{Float64},
    size_classes::Matrix{SizeClass},
    coral_spec::CoralSpec,
)
    plot_size_class(size_classes, t)
    k_area::Float64 = _k_area(cover[t-1, :, :], coral_spec.k_max)

    coral_prop = dropdims(sum(cover[t-1, :, :], dims=2), dims=2) ./ dropdims(sum(cover[t-1, :, :], dims=(1, 2)), dims=(1, 2))
    #_virtual_cover = virtual_cover(size_classes, coral_spec)
    settlers_cover::Vector{Float64} = coral_spec.settler_fracs .* coral_spec.k_max .* log(1 + k_area) .* coral_prop

    # Actual coral diameter growth [m/year]
    growth::Matrix{Float64} = linear_extension.(size_classes) .* log(1 + k_area)#k_area

    _, n_species, n_bins = size(cover)

    small = 1
    medium = collect(2:(n_bins-1))
    large = n_bins

    # Create new_size_classes
    # Small -> settlers
    small_cover_blocks = CoverBlock.(settlers_cover, interval_lower_bound.(size_classes[:, 1]), interval_upper_bound.(size_classes[:, 1]))
    small_size_classes = SizeClass.(
        small_cover_blocks,
        fill(small, n_species),
        linear_extension.(size_classes[:, 1]),
        survival_rate.(size_classes[:, 1])
    )

    medium_size_classes = new_medium_size_class.(
        size_classes[:, medium.-1],
        size_classes[:, medium],
        growth[:, medium.-1],
        growth[:, medium]
    )

    large_size_classes = new_large_size_class.(
        size_classes[:, large-1],
        size_classes[:, large],
        growth[:, large-1],
    )

    size_classes[:, small] = small_size_classes
    size_classes[:, medium] = medium_size_classes
    size_classes[:, large] = large_size_classes
    return size_class_cover.(size_classes)

    #@info "Cover after iteration $t-1: $(cover[t-1,:,:])"
    #@info "K-area $k"
end
function apply_changes!(size_class::Matrix{SizeClass}, reduction_density::SubArray{<:Union{Float32, Float64},2})::Nothing
    for i in axes(size_class, 1), j in axes(size_class, 2)
        apply_changes!.(size_class[i, j].cover_blocks, reduction_density[i, j])
    end
    return nothing
end
function apply_changes!(cover_block::CoverBlock, reduction_density::Union{Float32,Float64})::Nothing
    cover_block.diameter_density *= reduction_density

    return nothing
end

function timestep(
    cover::Array{Float64},
    recruits::Array{Float64},
    size_classes::Matrix{SizeClass},
    max_available::Float64,
    timestep::Int,
    plot::Bool=false
)
    if plot
        ext = Base.get_extension(parentmodule(blocks_model), :DynamicCoralPlotExt)
        if !isnothing(ext)
            parentmodule(blocks_model).Plot.plot_size_class(size_classes, timestep)
        end
    end
    k_area::Float64 = _k_area(cover, max_available)

    # Actual coral diameter growth [m/year]
    growth::Matrix{Float64} = linear_extension.(size_classes) .* log(1 + k_area)#k_area

    n_species, n_bins = size(cover)

    small = 1
    medium = collect(2:(n_bins-1))
    large = n_bins

    # Create new_size_classes
    # Small -> settlers
    small_size_classes = new_small_size_class.(
        size_classes[:, small],
        recruits,
        growth[:, small],
    )

    medium_size_classes = new_medium_size_class.(
        size_classes[:, medium.-1],
        size_classes[:, medium],
        growth[:, medium.-1],
        growth[:, medium]
    )

    large_size_classes = new_large_size_class.(
        size_classes[:, large-1],
        size_classes[:, large],
        growth[:, large-1],
    )

    size_classes[:, small] = small_size_classes
    size_classes[:, medium] = medium_size_classes
    size_classes[:, large] = large_size_classes
    return size_class_cover.(size_classes)

    #@info "Cover after iteration $t-1: $(cover[t-1,:,:])"
    #@info "K-area $k"
end
end
