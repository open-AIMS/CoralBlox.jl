using Test
using CoralBlox:
    max_projected_cover,
    _size_class_densities,
    _projected_cover,
    _adjusted_projected_cover,
    linear_extension_scale_factors

"""
Mock bin edges in m.
"""
function _mock_bin_edges(n_bins::Int64, n_functional_groups::Int64)::Matrix{Float64}
    return repeat(collect(range(0, 1, n_bins))', n_functional_groups, 1)
end

"""
Mock linear extensions in m.
"""
function _mock_linear_extensions(bin_edges::Matrix{Float64})::Matrix{Float64}
    return (bin_edges[:, 2:end] .- bin_edges[:, 1:(end - 1)]) .* 0.8
end

function mock_C_cover_t(
    _habitable_areas::Array{Float64, 3},
    n_functional_groups::Int64,
    n_size_classes::Int64,
    n_locs::Int64,
)
    return _habitable_areas .* reshape(
        hcat(
            [
                fill(i, n_functional_groups, n_size_classes) for i ∈ range(0.1, 0.9, n_locs)
            ]...,
        ),
        n_functional_groups,
        n_size_classes,
        n_locs,
    )
end

@testset "Adjutesd linear extension" begin
    n_functional_groups = 5
    n_size_classes = 7
    n_bins = n_size_classes + 1
    n_locs = 5

    _bin_edges::Matrix{Float64} = _mock_bin_edges(n_bins, n_functional_groups)
    _linear_extensions::Matrix{Float64} = _mock_linear_extensions(_bin_edges)
    _habitable_areas::Array{Float64, 3} =
        rand(n_functional_groups, n_size_classes, n_locs) .* 1e5
    _loc_habitable_areas::Vector{Float64} = dropdims(
        sum(_habitable_areas; dims=(1, 2)); dims=(1, 2)
    )
    _C_cover_t::Array{Float64, 3} = mock_C_cover_t(
        _habitable_areas, n_functional_groups, n_size_classes, n_locs
    )

    max_proj_cover::Vector{Float64} = CoralBlox.max_projected_cover(
        _linear_extensions, _bin_edges, _loc_habitable_areas
    )

    @testset "max_projected_cover" begin
        # Negative max_projected_cover values found
        @test all(max_proj_cover .>= 0)
        @test sum(max_proj_cover .< _loc_habitable_areas) == 0
    end

    size_class_densities::Array{Float64, 3} = CoralBlox._size_class_densities(
        _C_cover_t, _bin_edges
    )
    projected_cover::Vector{Float64} = CoralBlox._projected_cover(
        size_class_densities, _linear_extensions, _bin_edges
    )

    loc_cover::Vector{Float64} = dropdims(sum(_C_cover_t; dims=(1, 2)); dims=(1, 2))
    adjusted_projected_cover::Vector{Float64} = _adjusted_projected_cover(
        loc_cover, projected_cover, max_proj_cover, _loc_habitable_areas
    )

    @testset "projected_cover and adjusted_projected_cover" begin
        # Adjusted projected C(t) smaller than C(t-1)
        @test sum(loc_cover .> adjusted_projected_cover) == 0

        # Projected C(t) smaller than adjusted projected C(t)
        @test sum(adjusted_projected_cover .> projected_cover) == 0
    end

    @testset "linear_extension_scale_factors" begin
        scale_factors::AbstractVector{Float64} = CoralBlox.linear_extension_scale_factors(
            _C_cover_t, _loc_habitable_areas, _linear_extensions, _bin_edges, max_proj_cover
        )

        # Scale factors should be <= 1 and Scale factors should be >= 0
        @test sum(scale_factors .>= 1) == 0
        @test sum(scale_factors .<= 0) == 0
    end

    @testset "linear_extension_scale_factors — fully matured population (a == 0)" begin
        # All cover concentrated in the terminal size class (excluded from the target/growing
        # size classes), so every target block is empty and the quadratic solve's `a`
        # coefficient is exactly 0. This previously produced a NaN (0/0) or Inf scale factor.
        matured_cover = zeros(n_functional_groups, n_size_classes, n_locs)
        matured_cover[:, end, :] .= _habitable_areas[:, end, :]
        habitable_area = dropdims(sum(_habitable_areas; dims=(1, 2)); dims=(1, 2))

        scale_factors = CoralBlox.linear_extension_scale_factors(
            matured_cover, habitable_area, _linear_extensions, _bin_edges, max_proj_cover
        )

        @test all(isfinite, scale_factors)
        @test all(==(0.0), scale_factors)
    end
end
