module circular

# using DataStructures: CircularBuffer
using CircularArrayBuffers

struct SizeClass
    lower_bound::Float64
    upper_bound::Float64

    block_attrs::CircularArrayBuffer{Float64}
end

block_lower_bound(sc::SizeClass, block_id::Int64)::Float64 = sc.block_attrs[1, block_id]
block_upper_bound(sc::SizeClass, block_id::Int64)::Float64 = sc.block_attrs[2, block_id]
block_density(sc::SizeClass, block_id::Int64)::Float64 = sc.block_attrs[3, block_id]

block_lower_bounds(sc::SizeClass)::Vector{Float64} = sc.block_attrs[1, :]
block_upper_bounds(sc::SizeClass)::Vector{Float64} = sc.block_attrs[2, :]
block_densities(sc::SizeClass)::Vector{Float64} = sc.block_attrs[3, :]

mutable struct TerminalClass
    lower_bound::Float64
    upper_bound::Float64
    density::Float64
end

struct FunctionalGroup
    size_classes::Vector{SizeClass}
    terminal_class::TerminalClass
end

"""
    reallocate!(cb::CircularBuffer, n)

Resize CircularBuffer to the maximum capacity of n elements.
If n is smaller than the current buffer length, the first n elements will be retained.

# References
1. https://juliacollections.github.io/DataStructures.jl/latest/circ_buffer/
"""
function reallocate!(cb::CircularArrayBuffer, n::Integer)
    if n != capacity(cb)
        buf_new = Array{eltype(cb)}(undef, 3, n)
        len_new = min(size(cb, 2), n)
        # for i in 1:len_new
        #     buf_new[:, i] = cb[:, i]
        # end
        buf_new[:, 1:len_new] .= cb[:, :]

        cb.first = 1
        cb.step_size = len_new
        cb.buffer = buf_new
    end

    return cb
end

function reuse_buffers!(
    functional_groups::Vector{FunctionalGroup},
    cover::Union{Matrix{Float64}, SubArray{Float64, 2}}
)::Vector{FunctionalGroup}
    reuse_buffers!.(functional_groups, eachrow(cover))
    return functional_groups
end
function reuse_buffers!(functional_group::FunctionalGroup, cover::AbstractVector{Float64})::FunctionalGroup
    reuse_buffers!.(functional_group.size_classes, cover[1:end-1])

    area_factor::Float64 = average_area(functional_group.terminal_class)
    density::Float64 = cover[end] / area_factor
    functional_group.terminal_class.density = density

    return functional_group
end
function reuse_buffers!(size_class::SizeClass, cover::Float64)::Nothing
    empty!(size_class.block_attrs)

    area_factor::Float64 = average_area(size_class)
    density::Float64 = cover / area_factor

    push!(size_class.block_attrs, (size_class.lower_bound, size_class.upper_bound, density))

    return nothing
end

function SizeClass(
    lower_bound::Float64,
    upper_bound::Float64,
    cover::Float64;
    capacity::Int64=10
)::SizeClass
    area_factor::Float64 = π / 12 * (upper_bound^3 - lower_bound^3)
    density::Float64 = cover / area_factor

    block_attrs::CircularArrayBuffer{Float64} = CircularArrayBuffer{Float64}(3, capacity)

    push!(block_attrs, (density, lower_bound, upper_bound))

    return SizeClass(
        lower_bound,
        upper_bound,
        block_attrs
    )
end

function CreateTerminalClass(
    lower_bound::Float64,
    upper_bound::Float64,
    cover::Float64
)::TerminalClass
    area_factor::Float64 = π / 12 * (upper_bound^3 - lower_bound^3)
    density::Float64 = cover / area_factor

    return TerminalClass(
        lower_bound,
        upper_bound,
        density
    )
end

function FunctionalGroup(
    lower_bounds::AbstractVector{Float64},
    upper_bounds::AbstractVector{Float64},
    cover::AbstractVector{Float64}
)::FunctionalGroup
    size_classes::Vector{SizeClass} = SizeClass.(lower_bounds[1:end-1], upper_bounds[1:end-1], cover[1:end-1])
    terminal_class::TerminalClass = CreateTerminalClass(lower_bounds[end], upper_bounds[end], cover[end])

    return FunctionalGroup(
        size_classes,
        terminal_class
    )
end

"""
    apply_survival!(functional_group::FunctionalGroup, survival_rate::Vector{Float64})::Nothing
    apply_survival!(size_class::SizeClass, survival_rate::Float64)::Nothing
    apply_survival!(functional_groups::Vector{FunctionalGroup}, survival_rate::Union{Matrix{Float64}, SubArray{Float64, 2}})::Nothing

Apply mortality/survival probability to coral densities in blocks.
"""
function apply_survival!(
    functional_groups::Vector{FunctionalGroup},
    survival_rate::Union{Matrix{Float64}, SubArray{Float64, 2}}
)::Nothing
    apply_survival!.(functional_groups, eachrow(survival_rate))
    return nothing
end
function apply_survival!(
    functional_group::FunctionalGroup,
    survival_rate::Union{Vector{Float64}, SubArray{Float64, 1}}
)::Nothing
    apply_survival!.(functional_group.size_classes, survival_rate[1:end-1])
    functional_group.terminal_class.density *= survival_rate[end]

    return nothing
end
function apply_survival!(size_class::SizeClass, survival_rate::Float64)::Nothing
    # size_class.block_densities.buffer .*= survival_rate

    # block_attrs follow pattern lower, upper, density, so 3 is density
    size_class.block_attrs[3, :] .*= survival_rate

    return nothing
end

"""
    _apply_internal_growth!(functional_group::FunctionalGroup, growth_rates::Vector{Float64})::Nothing
    _apply_internal_growth!(size_class::SizeClass, growth_rate::Float64)::Nothing

Move coral blocks within size classes. Constrain blocks within bounds of size class.
"""
function _apply_internal_growth!(
    functional_group::FunctionalGroup,
    growth_rates::Vector{Float64}
)::Nothing
    _apply_internal_growth!.(functional_group.size_classes, growth_rates)
    return nothing
end
function _apply_internal_growth!(size_class::SizeClass, growth_rate::Float64)::Nothing
    # Apply growth directly to underlying buffer
    size_class.block_attrs[1, :] .+= growth_rate

    # Grow upper bounds and clamp bounds by upper bound of size class
    class_upper_bound::Float64 = size_class.upper_bound
    for (block_idx, ub) in enumerate(block_upper_bounds(size_class))
        size_class.block_attrs[2, block_idx] = min(
            class_upper_bound, ub + growth_rate
        )
    end

    return nothing
end

"""
    n_blocks(size_class::SizeClass)::Int64

Calculate number of blocks in a size class.
"""
function n_blocks(size_class::SizeClass)::Int64
    return size(size_class.block_attrs, 2)
end

function n_corals(density::Float64, lower_bound::Float64, upper_bound::Float64)::Float64
    return density * (upper_bound - lower_bound)
end

"""
    added_block_density(size_class::SizeClass, terminal::TerminalClass, block_idx::Int64, growth_rate::Float64)::Float64

Calculate the density added to the terminal size class from the smaller size class.
"""
function added_block_density(
    size_class::SizeClass,
    terminal::TerminalClass,
    block_idx::Int64,
    growth_rate::Float64
)::Float64
    block_lb::Float64, block_ub::Float64, block_density::Float64 = size_class.block_attrs[:, block_idx]

    coral_count::Float64 = n_corals(
        block_density,
        max(block_lb + growth_rate, terminal.lower_bound),
        block_ub + growth_rate
    )
    added_density::Float64 = coral_count / (terminal.upper_bound - terminal.lower_bound)

    return added_density
end

"""
    remove_outgrown!(size_class::SizeClass)::Nothing

Remove cover blocks that have outgrown the size class.
"""
function remove_outgrown!(size_class::SizeClass)::Nothing
    upper_bound::Float64 = size_class.upper_bound
    n_to_remove::Int64 = 0
    for block_lb in block_lower_bounds(size_class)
        if block_lb < upper_bound
            break
        end
        n_to_remove += 1
    end

    for i in 1:n_to_remove
        popfirst!(size_class.block_attrs)
    end

    return nothing
end

"""
    adjusted_growth(bound::Float64, upper_bound::Float64, prev_growth_rate::Float64, next_growth_rate::Float64)::Float64

Calculate the change in a bound when crossing a size class boundary.
"""
function adjusted_growth(
    bound::Float64,
    upper_bound::Float64,
    prev_growth_rate::Float64,
    next_growth_rate::Float64
)::Float64
    time_to_bound::Float64 = (upper_bound - bound) / prev_growth_rate
    return time_to_bound * prev_growth_rate + (1 - time_to_bound) * next_growth_rate
end

"""
    add_block!(size_class::SizeClass, density::Float64)::Nothing
    add_block!(size_class::SizeClass, lower_bound::Float64, upper_bound::Float64, density::Float64)::Nothing

Add new cover block to the given size class. Resize the size class if the buffer is already
full.
"""
function add_block!(size_class::SizeClass, density::Float64)::Nothing
    add_block!(size_class, size_class.lower_bound, size_class.upper_bound, density)

    return nothing
end
function add_block!(
    size_class::SizeClass, lower_bound::Float64, upper_bound::Float64, density::Float64
)::Nothing
    # Reallocate excess memory if buffers are full
    if capacity(size_class.block_attrs) == n_blocks(size_class)
        reallocate!(size_class.block_attrs, capacity(size_class.block_attrs) + 100)
    end

    push!(size_class.block_attrs, (lower_bound, upper_bound, density))

    return nothing
end

"""
    calculate_new_block!(block_lb::Float64, block_ub::Float64, block_density::Float64, next_class::SizeClass, prev_growth_rate::Float64, next_growth_rate::Float64)::Nothing

Calculate the density, lower bound and upper bound of a block transfitioning from a smaller
size class to a larger size class.
"""
function calculate_new_block!(
    curr_attrs::SubArray{Float64},
    next_class::SizeClass,
    prev_growth_rate::Float64,
    next_growth_rate::Float64
)::Tuple{Float64, Float64, Float64}
    # Check if the lower bound outgrows the upper bound as well
    curr_lb, curr_ub, curr_density = curr_attrs
    outgrowing_lb::Bool = curr_lb > (next_class.lower_bound - prev_growth_rate)
    # Calculate bounds and density of new cover block
    new_lower_bound::Float64 = (
        !outgrowing_lb ? next_class.lower_bound : curr_lb + adjusted_growth(
            curr_lb, next_class.lower_bound, prev_growth_rate, next_growth_rate
        )
    )
    new_upper_bound::Float64 = curr_ub + adjusted_growth(
        curr_ub, next_class.lower_bound, prev_growth_rate, next_growth_rate
    )

    proportion_moving::Float64 = outgrowing_lb ? 1.0 : 1.0 - (
        (next_class.lower_bound - (curr_lb + prev_growth_rate)) / (curr_ub - curr_lb)
    )

    # New Density = (number of corals * proportion moving)
    n_corals_moving::Float64 = curr_density * (curr_ub - curr_lb) * proportion_moving
    new_density::Float64 = n_corals_moving / (new_upper_bound - new_lower_bound)

    return new_lower_bound, new_upper_bound, new_density
end

"""
    transfer_blocks!(prev_class::SizeClass, next_class::SizeClass, prev_growth_rate::Float64, next_growth_rate::Float64)::Nothing
    transfer_blocks!(prev_class::SizeClass, terminal::TerminalClass, growth_rate::Float64)::Nothing

Transfer blocks from smaller size class to larger size class/Terminal class.
"""
function transfer_blocks!(
    prev_class::SizeClass,
    next_class::SizeClass,
    prev_growth_rate::Float64,
    next_growth_rate::Float64
)::Nothing
    # Blocks that exceed this bound will move to the next size class
    moving_bound::Float64 = prev_class.upper_bound - prev_growth_rate

    for block_idx in 1:n_blocks(prev_class)
        # Skip blocks that are not migrating
        if block_upper_bound(prev_class, block_idx) <= moving_bound
            continue
        end

        new_lower_bound, new_upper_bound, new_density = calculate_new_block!(
            @view(prev_class.block_attrs[:, block_idx]),
            next_class,
            prev_growth_rate,
            next_growth_rate
        )
        add_block!(next_class, new_lower_bound, new_upper_bound, new_density)
    end

    return nothing
end
function transfer_blocks!(
    prev_class::SizeClass,
    terminal::TerminalClass,
    growth_rate::Float64
)::Nothing
    # Blocks that exceed this bound will move to the next size class
    moving_bound::Float64 = prev_class.upper_bound - growth_rate

    # Accumulate density to add to terminal class
    additional_density::Float64 = 0.0
    for block_idx in 1:n_blocks(prev_class)
        # Skip blocks that are not migrating
        if block_upper_bound(prev_class, block_idx) <= moving_bound
            continue
        end
        additional_density += added_block_density(
            prev_class, terminal, block_idx, growth_rate
        )
    end

    terminal.density += additional_density

    return nothing
end

function transfer_and_grow!(
    prev_class::SizeClass,
    next_class::SizeClass,
    prev_growth_rate::Float64,
    next_growth_rate::Float64,
)::Nothing

    transfer_blocks!(
        prev_class,
        next_class,
        prev_growth_rate,
        next_growth_rate
    )

    # Grow corals in size classes
    _apply_internal_growth!(prev_class, prev_growth_rate)

    # Remove corals that out grow bounds
    remove_outgrown!(prev_class)

    return nothing
end
function transfer_and_grow!(
    prev_class::SizeClass,
    terminal::TerminalClass,
    growth_rate::Float64,
)::Nothing

    transfer_blocks!(
        prev_class,
        terminal,
        growth_rate
    )

    # Grow corals in size classes
    _apply_internal_growth!(prev_class, growth_rate)

    # Remove corals that out grow bounds
    remove_outgrown!(prev_class)

    return nothing
end

function average_area(terminal::TerminalClass)::Float64
    return π / 12 * (terminal.upper_bound^3 - terminal.lower_bound^3)
end
function average_area(size_class::SizeClass)::Float64
    return π / 12 * (size_class.upper_bound^3 - size_class.lower_bound^3)
end

function timestep!(
    functional_group::FunctionalGroup,
    recruitment::Float64,
    growth_rate::Union{Vector{Float64}, SubArray{Float64, 1}},
    survival_rate::Union{Vector{Float64}, SubArray{Float64, 1}}
)::Nothing
    apply_survival!(functional_group, survival_rate)

    transfer_and_grow!(
        functional_group.size_classes[end],
        functional_group.terminal_class,
        growth_rate[end]
    )

    n_classes::Int64 = length(functional_group.size_classes)
    for size_idx in n_classes:-1:2
        transfer_and_grow!(
            functional_group.size_classes[size_idx-1],
            functional_group.size_classes[size_idx],
            growth_rate[size_idx-1],
            growth_rate[size_idx],
        )
    end

    area_factor::Float64 = average_area(functional_group.size_classes[1])
    density::Float64 = recruitment / area_factor
    add_block!(functional_group.size_classes[1], density)

    return nothing
end

function timestep!(
    functional_groups::Vector{FunctionalGroup},
    recruitment::Vector{Float64},
    growth_rate::Matrix{Float64},
    survival_rate::Matrix{Float64}
)::Nothing
    timestep!.(functional_groups, recruitment, eachrow(growth_rate), eachrow(survival_rate))

    return nothing
end

"""
    coral_cover(functional_group::Vector{FunctionalGroup}, cache::SubArray{Float64, 2})::Nothing
    coral_cover(functional_group::FunctionalGroup, cache::SubArray{Float64, 1})::Nothing
    coral_cover(size_class::SizeClass)::Float64
"""
function coral_cover(
    functional_group::Vector{FunctionalGroup},
    cache::SubArray{Float64, 2}
)::Nothing
    coral_cover.(functional_group, eachrow(cache))

    return nothing
end
function coral_cover(
    functional_group::FunctionalGroup,
    cache::SubArray{Float64, 1}
)::Nothing
    cache[1:end-1] .= coral_cover.(functional_group.size_classes)
    cache[end] = functional_group.terminal_class.density * average_area(
        functional_group.terminal_class
    )
    return nothing
end
function coral_cover(size_class::SizeClass)::Float64
    cover::Float64 = 0.0
    for i in 1:n_blocks(size_class)
        cover += block_density(size_class, i) * π / 12 * (
            block_upper_bound(size_class, i)^3 - block_lower_bound(size_class, i)^3
        )
    end

    return cover
end

end
