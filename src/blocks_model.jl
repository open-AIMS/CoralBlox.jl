module blocks_model

include("coral_spec.jl")
# TODO Use only Tuples? Maybe have an intermediate struct called CoralBlock or smth?
mutable struct CoverBlock{N<:NTuple{2,Float64}}
    diameter_density::Float64       # Number of corals / m
    interval::N
end

function CoverBlock(initial_cover::Float64, bin_i::Float64, bin_f::Float64)::CoverBlock
    return CoverBlock(initial_cover / area_factor(bin_i, bin_f), (bin_i, bin_f))
end

#TODO I should remove linear extension, size_class and survival_rate from here to improve performance
mutable struct SizeClass{N<:NTuple{2,Float64},C<:CoverBlock}
    cover_blocks::Vector{C}
    const interval::N
    const size_class::Int64
    const linear_extension::Float64           # m/year
    const survival_rate::Float64              # %
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
    return SizeClass(CoverBlock[cover_block], cover_block.interval, size_class, linear_extension, survival_rate)
end

function SizeClass(
    cover_blocks::Vector{CoverBlock},
    size_class::SizeClass
)
    return SizeClass(cover_blocks, size_class.interval, size_class.size_class, size_class.linear_extension, size_class.survival_rate)
end

""" Area factor for a corals with diameter between bin_i and bin_f """
area_factor(bin_i, bin_f)::Float64 = (π / 12.0) * (bin_f^3 - bin_i^3)

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
    @assert (sum(coral_spec.initial_cover_fracs) <= 1.0) "Sum of all relative cover should be ≤ 1"
    return coral_spec.k_max .* coral_spec.initial_cover_fracs
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
    new_cover_blocks = append!(CoverBlock[CoverBlock(recruits, current_size_class.interval[1], current_size_class.interval[2])], move_current_blocks(current_size_class, current_growth))
    return SizeClass(new_cover_blocks, current_size_class)
end
function new_small_size_class(
    current_size_class::SizeClass,
    current_growth::Float64
)::SizeClass
    new_cover_blocks = move_current_blocks!(current_size_class, current_growth)
    return SizeClass(new_cover_blocks, current_size_class)
end

function new_small_size_class!(
    current_size_class::SizeClass,
    recruits::Float64,
    current_growth::Float64
)::Nothing
    new_cover_blocks = append!(CoverBlock[CoverBlock(recruits, current_size_class.interval[1], current_size_class.interval[2])], move_current_blocks(current_size_class, current_growth))
    current_size_class.cover_blocks = new_cover_blocks

    return nothing
end
function new_small_size_class!(
    current_size_class::SizeClass,
    current_growth::Float64
)::Nothing
    current_size_class.cover_blocks = move_current_blocks!(current_size_class, current_growth)
    return nothing
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
function new_medium_size_class!(
    prev_size_class::SizeClass,
    current_size_class::SizeClass,
    prev_growth::Float64,
    current_growth::Float64
)::Nothing
    prev_blocks = move_prev_blocks!(prev_size_class, prev_growth)
    current_blocks = move_current_blocks!(current_size_class, current_growth)

    current_size_class.cover_blocks = vcat(prev_blocks, current_blocks)

    return nothing
end

function move_prev_blocks(size_class::SizeClass, growth::Float64)::Vector{CoverBlock}
    # Select blocks that are going to next SizeClass
    new_upper_bounds::Vector{Float64} = interval_upper_bound.(size_class.cover_blocks) .+ growth
    sc_upper_bound::Float64 = size_class.interval[2]
    blocks_within_range::BitVector = new_upper_bounds .> sc_upper_bound

    new_blocks = Vector{CoverBlock}(undef, sum(blocks_within_range))
    for (idx, cover_block) in enumerate(size_class.cover_blocks[blocks_within_range])
        grown_lb::Float64 = interval(cover_block)[1] + growth
        new_lower_bound::Float64 = grown_lb < sc_upper_bound ? sc_upper_bound : grown_lb
        new_upper_bound::Float64 = interval(cover_block)[2] + growth
        new_interval::Tuple{Float64, Float64} = (new_lower_bound, new_upper_bound)
        new_diameter_density::Float64 = cover_block.diameter_density * size_class.survival_rate
        new_blocks[idx] = CoverBlock(new_diameter_density, new_interval)
    end

    return new_blocks
end
function move_prev_blocks!(size_class::SizeClass, growth::Float64)::Vector{CoverBlock}
    # Select blocks that are going to next SizeClass
    new_upper_bounds::Vector{Float64} = interval_upper_bound.(size_class.cover_blocks) .+ growth
    sc_upper_bound::Float64 = size_class.interval[2]
    blocks_within_range::BitVector = new_upper_bounds .> sc_upper_bound

    for (idx, cover_block) in zip(findall(blocks_within_range), size_class.cover_blocks[blocks_within_range])
        grown_lb::Float64 = interval(cover_block)[1] + growth
        new_lower_bound::Float64 = grown_lb < sc_upper_bound ? sc_upper_bound : grown_lb
        new_upper_bound::Float64 = interval(cover_block)[2] + growth
        new_interval::Tuple{Float64, Float64} = (new_lower_bound, new_upper_bound)
        new_diameter_density::Float64 = cover_block.diameter_density * size_class.survival_rate

        size_class.cover_blocks[idx].diameter_density = new_diameter_density
        size_class.cover_blocks[idx].interval = new_interval
    end

    return size_class.cover_blocks[blocks_within_range]
end

function move_current_blocks(size_class::SizeClass, growth::Float64)::Vector{CoverBlock}
    new_lower_bounds = interval_lower_bound.(size_class.cover_blocks) .+ growth
    sc_upper_bound = size_class.interval[2]
    blocks_within_range::BitVector = new_lower_bounds .< sc_upper_bound

    new_blocks = Vector{CoverBlock}(undef, sum(blocks_within_range))
    for (idx, cover_block) in enumerate(size_class.cover_blocks[blocks_within_range])
        grown_ub::Float64 = interval(cover_block)[2] + growth
        new_upper_bound::Float64 = grown_ub < sc_upper_bound ? grown_ub : sc_upper_bound
        new_lower_bound::Float64 = interval_lower_bound(cover_block) + growth
        new_interval::Tuple{Float64,Float64} = (new_lower_bound, new_upper_bound)
        new_diameter_density::Float64 = cover_block.diameter_density * size_class.survival_rate
        new_blocks[idx] = CoverBlock(new_diameter_density, new_interval)
    end

    return new_blocks
end
function move_current_blocks!(size_class::SizeClass, growth::Float64)::Vector{CoverBlock}
    new_lower_bounds = interval_lower_bound.(size_class.cover_blocks) .+ growth
    sc_upper_bound = size_class.interval[2]
    blocks_within_range::BitVector = new_lower_bounds .< sc_upper_bound

    for (idx, cover_block) in zip(findall(blocks_within_range), size_class.cover_blocks[blocks_within_range])
        grown_ub::Float64 = interval(cover_block)[2] + growth
        new_upper_bound::Float64 = grown_ub < sc_upper_bound ? grown_ub : sc_upper_bound
        new_lower_bound::Float64 = interval_lower_bound(cover_block) + growth
        new_interval::Tuple{Float64,Float64} = (new_lower_bound, new_upper_bound)
        new_diameter_density::Float64 = cover_block.diameter_density * size_class.survival_rate

        size_class.cover_blocks[idx].diameter_density = new_diameter_density
        size_class.cover_blocks[idx].interval = new_interval
    end

    return size_class.cover_blocks[blocks_within_range]
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
        new_diameter_density::Float64 = cover_block.diameter_density * prev_size_class.survival_rate
        interval_width::Float64 = cover_block.interval[2] + prev_growth - max(
            cover_block.interval[1] + prev_growth, current_size_class.interval[1]
        )
        n_new_corals += new_diameter_density * interval_width
    end

    current_Δinterval = current_size_class.interval[2] - current_size_class.interval[1]
    current_cover_block = current_size_class.cover_blocks[1]
    new_diameter_density = (n_new_corals / current_Δinterval) + current_cover_block.diameter_density * current_size_class.survival_rate

    new_cover_block = CoverBlock(new_diameter_density, current_cover_block.interval)
    return SizeClass(CoverBlock[new_cover_block], current_size_class)
end
function new_large_size_class!(
    prev_size_class::SizeClass,
    current_size_class::SizeClass,
    prev_growth::Float64,
)::Nothing
    new_upper_bound = interval_upper_bound.(prev_size_class.cover_blocks) .+ prev_growth
    blocks_within_range = new_upper_bound .> prev_size_class.interval[2]

    n_new_corals::Float64 = 0.0
    for cover_block in prev_size_class.cover_blocks[blocks_within_range]
        new_diameter_density::Float64 = cover_block.diameter_density * prev_size_class.survival_rate
        interval_width::Float64 = cover_block.interval[2] + prev_growth - max(
            cover_block.interval[1] + prev_growth, current_size_class.interval[1]
        )
        n_new_corals += new_diameter_density * interval_width
    end

    current_Δinterval = current_size_class.interval[2] - current_size_class.interval[1]
    current_cover_block = current_size_class.cover_blocks[1]
    new_diameter_density = (n_new_corals / current_Δinterval) + current_cover_block.diameter_density * current_size_class.survival_rate

    current_size_class.cover_blocks = CoverBlock[CoverBlock(new_diameter_density, current_cover_block.interval)]

    # new_cover_block = CoverBlock(new_diameter_density, current_cover_block.interval)
    # return SizeClass(CoverBlock[new_cover_block], current_size_class)

    return nothing
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

function apply_changes!(size_class::Matrix{SizeClass}, reduction_density::SubArray{<:Union{Float32,Float64},2})::Nothing
    for i in axes(size_class, 1), j in axes(size_class, 2)
        apply_changes!.(size_class[i, j].cover_blocks, reduction_density[i, j])
    end
    return nothing
end
function apply_changes!(cover_block::CoverBlock, reduction_density::Union{Float32,Float64})::Nothing
    cover_block.diameter_density *= reduction_density

    return nothing
end

"""
    timestep(cover::Matrix{Float64}, recruits::Vector{Float64}, size_classes::Matrix{SizeClass}, max_available_area::Float64, timestep::Int, plot::Bool=false)

Run one timestep iteration.

# Arguments
- `cover` : Coral cover with shape [groups ⋅ sizes]
- `recruits` : Recruits cover for each group
- `size_classes` : SizeClass instances with shape [groups ⋅ sizes]
- `max_available_area` : Maximum available area for corals to live
- `timestep` : Iteration timestep
- `plot` : If true, plot size classes
"""
function timestep(
    cover::Matrix{Float64},
    recruits::Vector{Float64},
    size_classes::Matrix{SizeClass},
    max_available_area::Float64,
    timestep::Int,
    plot::Bool=false
)
    if plot
        ext = Base.get_extension(parentmodule(blocks_model), :DynamicCoralPlotExt)
        if !isnothing(ext)
            parentmodule(blocks_model).Plot.plot_size_class(size_classes, timestep)
        end
    end

    # Coral diameter growth
    k_area::Float64 = _k_area(cover, max_available_area)
    growth::Matrix{Float64} = linear_extension.(size_classes) .* log(1 + k_area)

    _, n_bins = size(cover)

    small = 1
    medium = collect(2:(n_bins-1))
    large = n_bins

    # Create new_size_classes
    # small_size_classes = new_small_size_class.(
    #     size_classes[:, small],
    #     recruits,
    #     growth[:, small],
    # )
    new_small_size_class!.(
        size_classes[:, small],
        recruits,
        growth[:, small],
    )

    # medium_size_classes = new_medium_size_class.(
    #     size_classes[:, medium.-1],
    #     size_classes[:, medium],
    #     growth[:, medium.-1],
    #     growth[:, medium]
    # )
    new_medium_size_class!.(
        size_classes[:, medium.-1],
        size_classes[:, medium],
        growth[:, medium.-1],
        growth[:, medium]
    )

    # large_size_classes = new_large_size_class.(
    #     size_classes[:, large-1],
    #     size_classes[:, large],
    #     growth[:, large-1],
    # )
    new_large_size_class!.(
        size_classes[:, large-1],
        size_classes[:, large],
        growth[:, large-1],
    )

    # size_classes[:, small] = small_size_classes
    # size_classes[:, medium] = medium_size_classes
    # size_classes[:, large] = large_size_classes
    return size_class_cover.(size_classes)
end

end
