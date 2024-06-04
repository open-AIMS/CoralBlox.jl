using Test
using DynamicCoralCoverModel.circular:
    SizeClass,
    adjusted_growth,
    calculate_new_block!,
    remove_outgrown!,
    _apply_internal_growth!,
    transfer_blocks!,
    add_block!

@testset "adjusted_growth" begin

    adj_growth::Float64 = 0.0

    # Equal growth rates should return same growth rates
    growth_rate::Float64 = 2.0
    for _ in 1:50
        growth_rate = rand(1)[1]
        @test adjusted_growth(0.0, 1.0, growth_rate, growth_rate) ≈ growth_rate ||
            "Equal size class growth rates should yeild the same adjusted rate"
    end

    # If it takes 50% time to reach the boundary then growth is average of both rates
    growth_rate_prev::Float64 = 0.0
    growth_rate_next::Float64 = 0.0
    bound::Float64 = 0.0
    for _ in 1:50
        growth_rate_prev = rand(1)[1]
        growth_rate_next = rand(1)[1]
        lb = rand(1)[1]
        @test adjusted_growth(
            lb, lb + 0.5 * growth_rate_prev, growth_rate_prev, growth_rate_next
        ) ≈ (growth_rate_prev + growth_rate_next) * 0.5 ||
            "Crossing the boundary at time 0.5 should return average of both growth rates"
    end

    # Adjusted growth should be equal to next_growth if bound is equal to upper bound
    for _ in 1:50
        growth_rate_prev = rand(1)[1]
        growth_rate_next = rand(1)[1]
        bound = rand(1)[1]
        @test adjusted_growth(
            bound, bound, growth_rate_prev, growth_rate_next
        ) == growth_rate_next ||
            "If the bound is equal to the upper bound the growth must equal next growth"
    end

    # Adjusted growth rate should be a weighted average of both growth rates where the
    # weight is the time taken to cross the boundary
    proportion::Float64 = 0.0
    for _ in 1:50
        growth_rate_prev = rand(1)[1]
        growth_rate_next = rand(1)[1]
        lb = rand(1)[1]
        proportion = rand(1)[1]
        @test adjusted_growth(
            lb, lb + proportion * growth_rate_prev, growth_rate_prev, growth_rate_next
        ) ≈ growth_rate_prev * proportion + growth_rate_next * (1 - proportion) ||
            "Adjusted growth should return weighted edaverage of both growth rates"
    end
end

function equal_coral_counts(
    block_1::Tuple{Float64, Float64, Float64}, # lower, upper, density
    block_2::Tuple{Float64, Float64, Float64}
)::Bool
    count1::Float64 = (block_1[2] - block_1[1]) * block_1[3]
    count2::Float64 = (block_2[2] - block_2[1]) * block_2[3]
    return count1 ≈ count2
end

function proportion_of_coral_counts(
    block_1::Tuple{Float64, Float64, Float64}, # lower, upper, density
    block_2::Tuple{Float64, Float64, Float64},
    proportion::Float64
)::Bool
    count1::Float64 = (block_1[2] - block_1[1]) * block_1[3]
    count2::Float64 = (block_2[2] - block_2[1]) * block_2[3]
    return proportion * count1 ≈ count2
end

@testset "calculate_new_blocks" begin
    # The following tests assume adjusted_growth is correctly implemented

    next_class::SizeClass = SizeClass(1.0, 2.0, 1.0)
    # Coral counts shouid be preserved if the block completely leaves the lower size class
    for _ in 1:50
        class_lb = 1.0 + rand(1)[1]
        next_class = SizeClass(class_lb, 10.0, 5.0)

        block_lb = rand(1)[1]
        block_ub = class_lb
        block_density = 4.0 + rand(1)[1]

        growth_rate_prev = 2.0 + rand(1)[1]
        growth_rate_next = 2.0 + rand(1)[1]


        block_desc = calculate_new_block!(block_lb, block_ub, block_density, next_class, growth_rate_prev, growth_rate_next)
        @test equal_coral_counts((block_lb, block_ub, block_density), block_desc)
        @test block_desc[1] == block_lb + adjusted_growth(
            block_lb, next_class.lower_bound, growth_rate_prev, growth_rate_next
        )
        @test block_desc[2] == block_ub + adjusted_growth(
            block_ub, next_class.lower_bound, growth_rate_prev, growth_rate_next
        )
    end

    # Corals that migrate partially over the class bound should split the coral counts by
    # the proportion of the block transfered
    for _ in 1:50
        class_lb = 1.0 + rand(1)[1]
        next_class = SizeClass(class_lb, 10.0, 5.0)

        block_lb = rand(1)[1]
        block_ub = class_lb
        proportion = rand(1)[1]
        block_density = 4.0 * rand(1)[1]

        growth_rate_prev = (block_ub - block_lb) * proportion
        growth_rate_next = rand(1)[1]
        block_desc = calculate_new_block!(
            block_lb,
            block_ub,
            block_density,
            next_class,
            growth_rate_prev,
            growth_rate_next
        )
        @test proportion_of_coral_counts(
            (block_lb, block_ub, block_density),
            block_desc,
            proportion
        )
        @test block_desc[1] == class_lb
        @test block_desc[2] == block_ub + adjusted_growth(
            block_ub, next_class.lower_bound, growth_rate_prev, growth_rate_next
        )
    end
end

@testset "add_block!" begin
    test_size_class::SizeClass = SizeClass(2.0, 5.0, 10.0; capacity=5)

    lower_bound::Float64 = 2.0
    upper_bound::Float64 = 3.0
    density::Float64 = 1.0
    for i in 2:5
        add_block!(test_size_class, lower_bound, upper_bound, density)
        @test test_size_class.block_lower_bounds[i] == lower_bound ||
            "Block lower bounds not added to size class correctly"
        @test test_size_class.block_upper_bounds[i] == upper_bound ||
            "Block upper bounds not added to size class correctly"
        @test test_size_class.block_densities[i] == density ||
            "Block density not added to size_class correctly"
        @test (
            test_size_class.block_lower_bounds.length == i &&
            test_size_class.block_upper_bounds.length == i &&
            test_size_class.block_densities.length == i
        ) || "Size class should be of length $(i)"
    end

    @test (
        test_size_class.block_lower_bounds.capacity == 5 &&
        test_size_class.block_lower_bounds.capacity == 5 &&
        test_size_class.block_lower_bounds.capacity == 5
    ) || "Size class should still have capacity 5"

    # Trigger reallocation and check capacity
    add_block!(test_size_class, lower_bound, upper_bound, density)

    @test (
        test_size_class.block_lower_bounds.capacity == 15 &&
        test_size_class.block_lower_bounds.capacity == 15 &&
        test_size_class.block_lower_bounds.capacity == 15
    ) || "Size class should have capacity 15"

    # Test adding block that covers size class
    new_density::Float64 = 1.23
    add_block!(test_size_class, new_density)

    @test test_size_class.block_lower_bounds[end] == test_size_class.lower_bound ||
        "Block lower bound should be equal to size class lower bound"
    @test test_size_class.block_upper_bounds[end] == test_size_class.upper_bound ||
        "Block upper bound should be equal to size class upper bound"
    @test test_size_class.block_densities[end] == new_density ||
        "Block density should be equal to newest density"
end

function remove_first_block!(size_class::SizeClass)::Nothing
    if size_class.block_densities.length == 0
        throw(ArgumentError("Size Class does not contain any blocks to remove"))
    end
    popfirst!(size_class.block_lower_bounds)
    popfirst!(size_class.block_upper_bounds)
    popfirst!(size_class.block_densities)

    return nothing
end

@testset "transfer_blocks!" begin
    # Test equal growth rates
    prev_class::SizeClass = SizeClass(2.0, 5.0, 10.0)
    next_class::SizeClass = SizeClass(5.0, 10.0, 10.0)

    prev_growth_rate::Float64 = 1.0
    next_growth_rate::Float64 = 1.0

    transfer_blocks!(prev_class, next_class, prev_growth_rate, next_growth_rate)

    @test (
        next_class.block_densities.length == 2 &&
        next_class.block_lower_bounds.length == 2 &&
        next_class.block_upper_bounds.length == 2
    ) || "Length of Circuler Buffers should be 2 and in in sync"

    @test next_class.block_lower_bounds[2] == next_class.lower_bound ||
        "Lower bound of new block should be equal to lower_bound of class"
    @test next_class.block_upper_bounds[2] == next_class.lower_bound + next_growth_rate ||
        "Upper bound of new block should be equal to lower bound of class plus growth"
    # Equal size class growth rates gives equal density
    @test next_class.block_densities[2] ≈ prev_class.block_densities[1] ||
        "Block densities of should be equal given equal growth rates"

    # Test different growth rates
    prev_class = SizeClass(4.0, 5.0, 10.0)
    next_class = SizeClass(5.0, 10.0, 10.0)

    prev_growth_rate = 2.0
    next_growth_rate = 3.0

    expected_lower_bound::Float64 = 6.5
    expected_upper_bound::Float64 = 8.0

    transfer_blocks!(prev_class, next_class, prev_growth_rate, next_growth_rate)

    @test next_class.block_lower_bounds[2] == prev_class.lower_bound + 0.5 * (
        prev_growth_rate + next_growth_rate
    ) || "Lower bound of new block should be equal to lower_bound of class"
    @test next_class.block_upper_bounds[2] == next_class.lower_bound + next_growth_rate ||
        "Upper bound of new block should be equal to lower bound of class plus growth"

    # Test no transfers
    prev_class = SizeClass(3.0, 5.0, 10.0)
    next_class = SizeClass(5.0, 10.0, 10.0)

    remove_first_block!(prev_class)
    add_block!(prev_class, 3.0, 4.0, 1.23)

    # Assure previous growth rate is lower then the gap between block and next size class
    prev_growth_rate = 1.0
    next_growth_rate = 3.0

    transfer_blocks!(prev_class, next_class, prev_growth_rate, next_growth_rate)

    @test (
        prev_class.block_lower_bounds.length == 1 &&
        prev_class.block_upper_bounds.length == 1 &&
        prev_class.block_densities.length == 1
    ) || "Previous size class should only contain one block"

    @test (
        next_class.block_lower_bounds.length == 1 &&
        next_class.block_upper_bounds.length == 1 &&
        next_class.block_densities.length == 1
    ) || "Next size class should only contain one block. None should be added."
end

function equal_bounds(
    size_class::SizeClass,
    lower_bounds::Vector{Float64},
    upper_bounds::Vector{Float64}
)::Bool
    equal::Bool = true
    for (expected, lb) in zip(size_class.block_lower_bounds, lower_bounds)
        equal &= expected == lb
    end
    for (expected, ub) in zip(size_class.block_lower_bounds, upper_bounds)
        equal &= expected == ub
    end
    return equal
end

@testset "_apply_internal_growth!" begin
    test_size_class::SizeClass = SizeClass(1.0, 6.0, 1.0)
    remove_first_block!(test_size_class)

    lower_bounds::Vector{Float64} = reverse([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    upper_bounds::Vector{Float64} = reverse([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5])

    for i in 1:length(lower_bounds)
        add_block!(test_size_class, lower_bounds[i], upper_bounds[i], 2.0)
    end

    growth_rate::Float64 = 2.0
    expected_lower_bounds::Vector{Float64} = lower_bounds .+ growth_rate
    expected_upper_bounds::Vector{Float64} = reverse(
        [3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.0, 6.0, 6.0]
    )

    _apply_internal_growth!(test_size_class, growth_rate)

    @test all(test_size_class.block_lower_bounds .== expected_lower_bounds) ||
        "Unexpected block lower bounds."
    @test all(test_size_class.block_upper_bounds .== expected_upper_bounds) ||
        "Unexpected block upper bounds."
    @test all(test_size_class.block_densities .== 2.0) ||
        "Unexpected block densities."
end

@testset "remove_outgrown!" begin
    test_size_class::SizeClass = SizeClass(1.0, 6.0, 1.0)
    remove_first_block!(test_size_class)

    lower_bounds::Vector{Float64} = reverse([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    upper_bounds::Vector{Float64} = reverse([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5])

    for i in 1:length(lower_bounds)
        add_block!(test_size_class, lower_bounds[i], upper_bounds[i], 2.0)
    end

    growth_rate::Float64 = 2.0
    expected_lower_bounds::Vector{Float64} = (lower_bounds .+ growth_rate)[4:end]
    expected_upper_bounds::Vector{Float64} = reverse([3.5, 4.0, 4.5, 5.0, 5.5, 6.0])

    _apply_internal_growth!(test_size_class, growth_rate)
    remove_outgrown!(test_size_class)

    @test all(test_size_class.block_lower_bounds .== expected_lower_bounds) ||
        "Unexpected block lower bounds."
    @test all(test_size_class.block_upper_bounds .== expected_upper_bounds) ||
        "Unexpected block upper bounds."
    @test all(test_size_class.block_densities .== 2.0) ||
        "Unexpected block densities."
end
