using DataStructures: CircularBuffer

struct SizeClass
    lower_bound::Float64
    upper_bound::Float64

    block_lower_bounds::CircularBuffer{Float64}
    block_upper_bounds::CircularBuffer{Float64}
    block_densities::CircularBuffer{Float64}

    # Caches
    buf_new::Vector{Float64}
    movement_cache::Vector{Float64}  # 3 elements
end


function SizeClass(
    lower_bound::Float64,
    upper_bound::Float64,
    cover::Float64;
    capacity::Int64=256
)::SizeClass
    area_factor::Float64 = π / 12 * (upper_bound^3 - lower_bound^3)
    density::Float64 = cover / area_factor

    block_lower_bounds::CircularBuffer{Float64} = CircularBuffer{Float64}(capacity)
    block_upper_bounds::CircularBuffer{Float64} = CircularBuffer{Float64}(capacity)
    block_densities::CircularBuffer{Float64} = CircularBuffer{Float64}(capacity)

    # Assert statement will not be compiled/executed on higher optimisation levels
    @assert (
        length(block_densities) == length(block_lower_bounds) == length(block_upper_bounds)
    ) "Length of bounds and densities are not the same."

    push!(block_lower_bounds, lower_bound)
    push!(block_upper_bounds, upper_bound)
    push!(block_densities, density)

    buf_new::Vector{Float64} = zeros(capacity)
    cache::Vector{Float64} = zeros(3)
    return SizeClass(
        lower_bound,
        upper_bound,
        block_lower_bounds,
        block_upper_bounds,
        block_densities,
        buf_new,
        cache
    )
end

function Base.show(io::IO, mime::MIME"text/plain", sc::SizeClass)::Nothing
    println("""
    Size Class

    Lower Bound: $(sc.lower_bound)
    Upper Bound: $(sc.upper_bound)
    Current buffer size: $(length(sc.buf_new))
    """)

    return nothing
end

mutable struct TerminalClass
    lower_bound::Float64
    upper_bound::Float64
    density::Float64

    function TerminalClass(lower_bound::Float64, upper_bound::Float64, cover::Float64)
        area_factor::Float64 = π / 12 * (upper_bound^3 - lower_bound^3)
        density::Float64 = cover / area_factor

        return new(lower_bound, upper_bound, density)
    end
end


struct FunctionalGroup
    size_classes::Vector{SizeClass}
    terminal_class::TerminalClass
end

function FunctionalGroup(
    lower_bounds::AbstractVector{Float64},
    upper_bounds::AbstractVector{Float64},
    cover::AbstractVector{Float64}
)::FunctionalGroup
    size_classes::Vector{SizeClass} = SizeClass.(lower_bounds[1:end-1], upper_bounds[1:end-1], cover[1:end-1])
    terminal_class::TerminalClass = TerminalClass(lower_bounds[end], upper_bounds[end], cover[end])

    return FunctionalGroup(
        size_classes,
        terminal_class
    )
end

function Base.show(io::IO, mime::MIME"text/plain", fg::FunctionalGroup)::Nothing
    lb = fg.terminal_class.lower_bound
    ub = fg.terminal_class.upper_bound
    density = fg.terminal_class.density

    println("""
    Functional Group

    Number of Size Classes: $(length(fg.size_classes))

    Terminal Size Class:
        Lower bound: $(lb)
        Upper bound: $(ub)
        Current density: $(density)
        Current number of corals: $(n_corals(lb, ub, density))
    """)

    return nothing
end

"""
    reallocate!(cb::CircularBuffer, buf_new::Vector{Float64})

Resize CircularBuffer to the maximum capacity of n elements.
If n is smaller than the current buffer length, the first n elements will be retained.

# References
1. DataStructures.jl
"""
function reallocate!(cb::CircularBuffer, buf_new::Vector{Float64})
    n = length(buf_new)
    if n != cb.capacity
        len_new = min(length(cb), n)
        buf_new[1:len_new] .= cb[1:len_new]

        cb.capacity = n
        cb.first = 1
        cb.length = len_new
        cb.buffer = copy(buf_new)
    end
    return cb
end

function reuse_buffers!(
    functional_groups::Vector{FunctionalGroup},
    cover::Union{Matrix{Float64},SubArray{Float64,2}}
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
    empty!(size_class.block_lower_bounds)
    empty!(size_class.block_upper_bounds)
    empty!(size_class.block_densities)

    area_factor::Float64 = average_area(size_class)
    density::Float64 = cover / area_factor

    push!(size_class.block_lower_bounds, size_class.lower_bound)
    push!(size_class.block_upper_bounds, size_class.upper_bound)
    push!(size_class.block_densities, density)

    return nothing
end

"""
    apply_mortality!(functional_group::FunctionalGroup, survival_rate::Vector{Float64})::Nothing
    apply_mortality!(size_class::SizeClass, survival_rate::Float64)::Nothing
    apply_mortality!(functional_groups::Vector{FunctionalGroup}, survival_rate::Union{Matrix{Float64}, SubArray{Float64, 2}})::Nothing

Apply mortality/survival probability to coral densities in blocks.
"""
function apply_mortality!(
    functional_groups::Vector{FunctionalGroup},
    survival_rate::Union{Matrix{Float64},SubArray{Float64,2}}
)::Nothing
    apply_mortality!.(functional_groups, eachrow(survival_rate))
    return nothing
end
function apply_mortality!(
    functional_group::FunctionalGroup,
    survival_rate::Union{Vector{Float64},SubArray{Float64,1}}
)::Nothing
    apply_mortality!.(functional_group.size_classes, survival_rate[1:end-1])
    functional_group.terminal_class.density *= survival_rate[end]

    return nothing
end
function apply_mortality!(size_class::SizeClass, survival_rate::Float64)::Nothing
    size_class.block_densities .*= survival_rate

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
    size_class.block_lower_bounds .+= growth_rate
    size_class.block_upper_bounds .+= growth_rate

    # Grow upper bounds and clamp bounds by upper bound of size class
    class_upper_bound::Float64 = size_class.upper_bound
    for block_idx in 1:n_blocks(size_class)
        if size_class.block_upper_bounds[block_idx] < class_upper_bound
            break
        end
        size_class.block_upper_bounds[block_idx] = class_upper_bound
    end

    return nothing
end

"""
    n_blocks(size_class::SizeClass)::Int64

Calculate number of blocks in a size class.
"""
function n_blocks(size_class::SizeClass)::Int64
    return length(size_class.block_upper_bounds)
end

function n_corals(lower_bound::Float64, upper_bound::Float64, density::Float64)::Float64
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
    block_lb::Float64 = size_class.block_lower_bounds[block_idx]
    block_ub::Float64 = size_class.block_upper_bounds[block_idx]
    block_density::Float64 = size_class.block_densities[block_idx]

    coral_count::Float64 = n_corals(
        max(block_lb + growth_rate, terminal.lower_bound),
        block_ub + growth_rate,
        block_density
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
    for block_lb in size_class.block_lower_bounds
        if block_lb < upper_bound
            break
        end
        n_to_remove += 1
    end

    for _ in 1:n_to_remove
        popfirst!(size_class.block_lower_bounds)
        popfirst!(size_class.block_upper_bounds)
        popfirst!(size_class.block_densities)
    end

    return nothing
end

"""
    crossedge_displacement(bound::Float64, upper_bound::Float64, prev_growth_rate::Float64, next_growth_rate::Float64)::Float64

Calculate the change in a bound when crossing a size class edge.
"""
function crossedge_displacement(
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
function add_block!(size_class::SizeClass, block_attrs::SubArray{Float64,1})
    # Reallocate excess memory if buffers are full
    if size_class.block_densities.capacity == size_class.block_densities.length
        new_current_capacity::Int64 = size_class.block_densities.capacity + 32
        if length(size_class.buf_new) != new_current_capacity
            resize!(size_class.buf_new, new_current_capacity)
        end

        reallocate!(size_class.block_lower_bounds, size_class.buf_new)
        reallocate!(size_class.block_upper_bounds, size_class.buf_new)
        reallocate!(size_class.block_densities, size_class.buf_new)
    end

    push!(size_class.block_lower_bounds, block_attrs[1])
    push!(size_class.block_upper_bounds, block_attrs[2])
    push!(size_class.block_densities, block_attrs[3])

    return nothing
end
function add_block!(size_class::SizeClass, density::Float64)::Nothing
    add_block!(size_class, size_class.lower_bound, size_class.upper_bound, density)

    return nothing
end
function add_block!(
    size_class::SizeClass, lower_bound::Float64, upper_bound::Float64, density::Float64
)::Nothing
    # Reallocate excess memory if buffers are full
    if size_class.block_densities.capacity == size_class.block_densities.length
        new_current_capacity::Int64 = size_class.block_densities.capacity + 32
        if length(size_class.buf_new) != new_current_capacity
            resize!(size_class.buf_new, new_current_capacity)
        end

        reallocate!(size_class.block_lower_bounds, size_class.buf_new)
        reallocate!(size_class.block_upper_bounds, size_class.buf_new)
        reallocate!(size_class.block_densities, size_class.buf_new)
    end

    push!(size_class.block_lower_bounds, lower_bound)
    push!(size_class.block_upper_bounds, upper_bound)
    push!(size_class.block_densities, density)

    return nothing
end

"""
    calculate_new_block!(block_lb::Float64, block_ub::Float64, block_density::Float64, next_class::SizeClass, prev_growth_rate::Float64, next_growth_rate::Float64)::Nothing

Calculate the density, lower bound and upper bound of a block transitioning from a smaller
size class to a larger size class.
"""
function calculate_new_block!(
    block_lb::Float64,
    block_ub::Float64,
    block_density::Float64,
    next_class::SizeClass,
    prev_growth_rate::Float64,
    next_growth_rate::Float64
)::Tuple{Float64,Float64,Float64}
    # Check if the lower bound outgrows the upper bound as well
    outgrowing_lb::Bool = block_lb > (next_class.lower_bound - prev_growth_rate)

    # Calculate bounds and density of new cover block
    new_lower_bound::Float64 = (
        !outgrowing_lb ? next_class.lower_bound : block_lb + crossedge_displacement(
            block_lb, next_class.lower_bound, prev_growth_rate, next_growth_rate
        )
    )

    new_upper_bound::Float64 = block_ub + crossedge_displacement(
        block_ub, next_class.lower_bound, prev_growth_rate, next_growth_rate
    )

    proportion_moving::Float64 = outgrowing_lb ? 1.0 : 1.0 - (
        (next_class.lower_bound - (block_lb + prev_growth_rate)) / (block_ub - block_lb)
    )

    # New Density = (number of corals * proportion moving)
    n_corals_moving::Float64 = block_density * (block_ub - block_lb) * proportion_moving
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
        if prev_class.block_upper_bounds[block_idx] <= moving_bound
            break
        end
        prev_class.movement_cache .= calculate_new_block!(
            prev_class.block_lower_bounds[block_idx],
            prev_class.block_upper_bounds[block_idx],
            prev_class.block_densities[block_idx],
            next_class,
            prev_growth_rate,
            next_growth_rate
        )
        add_block!(next_class, @view(prev_class.movement_cache[1:3]))
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
        if prev_class.block_upper_bounds[block_idx] <= moving_bound
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
    if prev_growth_rate == 0.0
        return nothing
    end

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
    if growth_rate == 0.0
        return nothing
    end

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

"""
    merge_transfer!(smallest_class::SizeClass, next_class::SizeClass, smallest_growth_rate::Float64, next_growth_rate::Float64)::Nothing

All cover blocks within the smallest class will always transfer part of the block to the
next size class.
"""
function merge_transfer!(
    smallest_class::SizeClass,
    next_class::SizeClass,
    smallest_growth_rate::Float64,
    next_growth_rate::Float64
)::Nothing
    lower_bound_condition::Float64 = smallest_class.upper_bound - smallest_growth_rate

    # Migrate blocks that do not give identical new blocks
    final_index::Int64 = 0
    early_exit::Bool = false
    for block_idx in 1:n_blocks(smallest_class)
        if smallest_class.block_lower_bounds[block_idx] <= lower_bound_condition
            final_index = block_idx
            early_exit = true
            break
        end
        smallest_class.movement_cache .= calculate_new_block!(
            smallest_class.block_lower_bounds[block_idx],
            smallest_class.block_upper_bounds[block_idx],
            smallest_class.block_densities[block_idx],
            next_class,
            smallest_growth_rate,
            next_growth_rate
        )
        add_block!(
            next_class,
            @view(smallest_class.movement_cache[1:3])
        )
    end

    if !early_exit
        return nothing
    end
    # The width of new blocks is always next_growth_rate. So we need only add densities and
    # adjust for new width
    new_density::Float64 = sum(smallest_class.block_densities[final_index:end])
    new_density *= smallest_growth_rate / next_growth_rate

    add_block!(
        next_class,
        smallest_class.upper_bound,
        smallest_class.upper_bound + next_growth_rate,
        new_density
    )
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
    recruits::Float64,
    growth_rate::SubArray{Float64,1},
    survival_rate::SubArray{Float64,1}
)::Nothing
    apply_mortality!(functional_group, survival_rate)

    transfer_and_grow!(
        functional_group.size_classes[end],
        functional_group.terminal_class,
        growth_rate[end]
    )

    n_classes::Int64 = length(functional_group.size_classes)
    for size_idx in n_classes:-1:3
        transfer_and_grow!(
            functional_group.size_classes[size_idx-1],
            functional_group.size_classes[size_idx],
            growth_rate[size_idx-1],
            growth_rate[size_idx],
        )
    end

    if growth_rate[1] != 0.0
        merge_transfer!(
            functional_group.size_classes[1],
            functional_group.size_classes[2],
            growth_rate[1],
            growth_rate[2]
        )
        _apply_internal_growth!(functional_group.size_classes[1], growth_rate[1])
        remove_outgrown!(functional_group.size_classes[1])
    end

    # Add recruits cover block only when there are recruits
    if recruits > 0.0
        area_factor::Float64 = average_area(functional_group.size_classes[1])
        recruits_density::Float64 = recruits / area_factor
        add_block!(functional_group.size_classes[1], recruits_density)
    end

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
    coral_cover(functional_group::Vector{FunctionalGroup}, C_cover::SubArray{Float64, 2})::Nothing
    coral_cover(functional_group::FunctionalGroup, C_cover::SubArray{Float64, 1})::Nothing
    coral_cover(size_class::SizeClass)::Float64
"""
function coral_cover(
    functional_group::Vector{FunctionalGroup},
    C_cover::SubArray{Float64,2}
)::Nothing
    coral_cover.(functional_group, eachrow(C_cover))

    return nothing
end
function coral_cover(
    functional_group::FunctionalGroup,
    C_cover::SubArray{Float64,1}
)::Nothing
    C_cover[1:end-1] .= coral_cover.(functional_group.size_classes)
    C_cover[end] = functional_group.terminal_class.density * average_area(
        functional_group.terminal_class
    )
    return nothing
end
function coral_cover(size_class::SizeClass)::Float64
    cover::Float64 = 0.0
    #? I think we could broadcast this operation
    for i in 1:n_blocks(size_class)
        cover += size_class.block_densities[i] * π / 12 * (
            size_class.block_upper_bounds[i]^3 - size_class.block_lower_bounds[i]^3
        )
    end

    return cover
end
