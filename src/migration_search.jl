using DataStructures: CircularBuffer

"""
    _n_migrating(block_upper_bounds::CircularBuffer{Float64}, moving_bound::Float64)::Int64

Number of leading blocks whose upper bound is greater than `moving_bound` — the count of
migrating blocks for `transfer_blocks!`.

# Background

A block migrates to the next size class iff, after growing by `growth_rate`, its upper
edge crosses the class ceiling:

    block_upper_bound + growth_rate > class.upper_bound
  ⟺ block_upper_bound              > class.upper_bound - growth_rate   (== moving_bound)

`block_upper_bounds` is kept sorted descending, so the migrating blocks are always the
prefix `1:n`; this returns `n`, the boundary between `> moving_bound` and `<= moving_bound`.

`isless(moving_bound, x)` — rather than `x > moving_bound` — is used so the result matches
`searchsortedfirst(block_upper_bounds, moving_bound; rev=true) - 1` exactly, including ties
and signed zeros.

# Why a linear scan

Earlier versions used `searchsortedfirst(...; rev=true)`. Profiling ADRIA showed that call
was ~35% of CoralBlox stepping time: the keyword form allocates a `ReverseOrdering` per
call and every probe goes through `CircularBuffer`'s modular `getindex`.

Four strategies were compared — `searchsortedfirst`, this linear scan, and two hand-rolled
binary searches over the backing `Vector` (see `benchmark/migration_search.jl`). A
micro-benchmark favoured binary search, but in the full model it and the linear scan came
out even: both cut the cutoff-search cost ~25-30% and `transfer_blocks!` ~15% versus
`searchsortedfirst` (~8% off total ADRIA runtime). Binary search wins only when many
blocks migrate; here `moving_bound` sits just below the class ceiling so `n` is typically
1-3, and a couple of sequential, prefetchable reads with an early break match it. The
linear scan was chosen for being simpler and obviously correct at equal performance.
"""
function _n_migrating(
    block_upper_bounds::CircularBuffer{Float64}, moving_bound::Float64
)::Int64
    n::Int64 = 0
    @inbounds for i ∈ 1:length(block_upper_bounds)
        isless(moving_bound, block_upper_bounds[i]) || break
        n += 1
    end
    return n
end
