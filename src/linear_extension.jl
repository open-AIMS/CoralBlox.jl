"""
    max_projected_cover(linear_extensions::Matrix{Float64}, bin_edges::Matrix{Float64}, habitable_area::Float64)::Float64
    max_projected_cover(linear_extensions::Matrix{Float64}, bin_edges::Matrix{Float64}, habitable_areas::Vector{Float64})

Finds the maximum projected cover for each location considering that all habitable area is
concentrated in the in the last size class of the "worst case scenario" functional group,
and that it grows by that functional group's linear extension.

# Arguments
- `linear_extensions` : Functional Groups x Size Classes
- `bin_edges` : Functional Groups x Bin Edges
- `habitable_areas` : Habitable area
"""
function max_projected_cover(
    linear_extensions::Matrix{Float64},
    bin_edges::Matrix{Float64},
    habitable_area::Float64
)::Float64
    return _diameter_coef(linear_extensions, bin_edges) * habitable_area
end
function max_projected_cover(
    linear_extensions::Matrix{Float64},
    bin_edges::Matrix{Float64},
    habitable_areas::Vector{Float64}
)::Vector{Float64}
    return [_diameter_coef(linear_extensions, bin_edges)] .* habitable_areas
end

function _diameter_coef(
    linear_extensions::Matrix{Float64}, bin_edges::Matrix{Float64}
)::Float64
    last_linear_extensions::Vector{Float64} = linear_extensions[:, end-1]
    last_C_bins::Vector{Float64} = bin_edges[:, end-1]
    max_diameter_idx::Int64 = findmax(last_linear_extensions .+ last_C_bins)[2]
    last_upper_bound::Float64 = last_C_bins[max_diameter_idx]
    last_linear_extension::Float64 = last_linear_extensions[max_diameter_idx]
    return ((last_upper_bound + last_linear_extension)^2) / (last_upper_bound^2)
end

"""
    adjusted_linear_extension(C_cover_t::AbstractArray{Float64,3}, loc_habitable_areas::AbstractVector{Float64}, linear_extensions::AbstractMatrix{Float64}, bin_edges::AbstractMatrix{Float64}, max_projected_cover::AbstractVector{Float64})

Adjusted linear extension.

# Arguments
- `C_cover_t` : Cover at previous timestep
- `loc_habitable_areas` : dimensions 1 x 1 x habitable_locations
- `linear_extensions` : dimensions functional_groups x size_classes
- `bin_edges` : dimensions functional_groups x bin_edges
- `max_projected_cover` : dimensions habitable_locations
"""
function linear_extension_scale_factors(
    C_cover_t::AbstractMatrix{Float64},
    habitable_area::Float64,
    linear_extensions::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64},
    max_projected_cover::Float64
)::Float64
    total_cover::Float64 = sum(C_cover_t)

    # Average density for each functional group and size class
    size_class_densities::Matrix{Float64} = _size_class_densities(C_cover_t, bin_edges)

    # Projected coral cover at t+1
    projected_cover::Float64 = _projected_cover(
        size_class_densities, linear_extensions, bin_edges
    )

    adjusted_projected_cover::Float64 = _adjusted_projected_cover(
        total_cover, projected_cover, max_projected_cover, habitable_area
    )

    # Solve quadratic equation
    a::Float64 = _quadratic_coeff(size_class_densities, linear_extensions, Δd(bin_edges, 1))
    b::Float64 = _linear_coeff(size_class_densities, linear_extensions, Δd(bin_edges, 2))
    c::Float64 = _constant_coeff(
        size_class_densities, adjusted_projected_cover, Δd(bin_edges, 3)
    )

    return (sqrt((b^2) - (4 * a * c)) - (b)) / (2 * a)
end
function linear_extension_scale_factors(
    C_cover_t::AbstractArray{Float64,3},
    loc_habitable_areas::AbstractVector{Float64},
    linear_extensions::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64},
    max_projected_cover::AbstractVector{Float64}
)::AbstractVector{Float64}
    # Total cover for each location
    loc_C_cover::Vector{Float64} = dropdims(sum(C_cover_t, dims=(1, 2)); dims=(1, 2))

    # Average density for each functional group, size class and location
    size_class_densities::Array{Float64,3} = _size_class_densities(C_cover_t, bin_edges)

    # Projected coral cover for each location at t+1
    projected_cover::Vector{Float64} = _projected_cover(
        size_class_densities, linear_extensions, bin_edges
    )

    adjusted_projected_cover::Vector{Float64} = _adjusted_projected_cover(
        loc_C_cover, projected_cover, max_projected_cover, loc_habitable_areas
    )

    # Solve quadratic equation
    a::Vector{Float64} = _quadratic_coeff(
        size_class_densities, linear_extensions, Δd(bin_edges, 1)
    )
    b::Vector{Float64} = _linear_coeff(
        size_class_densities, linear_extensions, Δd(bin_edges, 2)
    )
    c::Vector{Float64} = _constant_coeff(
        size_class_densities, adjusted_projected_cover, Δd(bin_edges, 3)
    )

    return (sqrt.((b .^ 2) .- (4 .* a .* c)) .- (b)) ./ (2 .* a)
end

Δd(
    bin_edges::AbstractMatrix{Float64},
    n::Int64
)::Matrix{Float64} = ((bin_edges[:, 2:end] .^ n) .- (bin_edges[:, 1:end-1] .^ n))

_size_class_densities(
    C_cover_t::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64}
)::Matrix{Float64} = (12 / π) .* (C_cover_t ./ Δd(bin_edges, 3))
_size_class_densities(
    C_cover_t::AbstractArray{Float64,3},
    bin_edges::AbstractMatrix{Float64}
)::Array{Float64,3} = (12 / π) .* (C_cover_t ./ Δd(bin_edges, 3))

function _projected_cover(
    size_class_densities::AbstractMatrix{Float64},
    linear_extensions::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64}
)::Float64
    return sum(
        (π / 12) .*
        size_class_densities .* (
            (3 .* (linear_extensions .^ 2) .* Δd(bin_edges, 1)) .+
            (3 .* linear_extensions .* Δd(bin_edges, 2)) .+
            Δd(bin_edges, 3)
        )
    )
end
function _projected_cover(
    size_class_densities::AbstractArray{Float64,3},
    linear_extensions::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64}
)::Vector{Float64}
    return dropdims(sum(
            ((π / 12) .*
             size_class_densities .* (
                (3 .* (linear_extensions .^ 2) .* Δd(bin_edges, 1)) .+
                (3 .* linear_extensions .* Δd(bin_edges, 2)) .+
                Δd(bin_edges, 3)
            )
            ),
            dims=(2, 1)
        ); dims=(2, 1))
end

function _adjusted_projected_cover(
    total_cover::Float64,
    projected_cover::Float64,
    max_projected_cover::Float64,
    habitable_area::Float64
)::Float64
    return total_cover + (
        (projected_cover - total_cover) *
        (
            (habitable_area - total_cover) / (max_projected_cover - total_cover)
        )
    )
end
function _adjusted_projected_cover(
    loc_cover::AbstractVector{Float64},
    projected_cover::AbstractVector{Float64},
    max_projected_cover::AbstractVector{Float64},
    loc_habitable_areas::AbstractVector{Float64}
)::Vector{Float64}
    return loc_cover .+ (
        (projected_cover .- loc_cover) .*
        (
            (loc_habitable_areas .- loc_cover) ./ (max_projected_cover .- loc_cover)
        )
    )
end

function _quadratic_coeff(
    size_class_densities::AbstractMatrix{Float64},
    linear_extensions::AbstractMatrix{Float64},
    Δ¹d::AbstractMatrix{Float64}
)::Float64
    return sum(3 .* (size_class_densities .* (linear_extensions .^ 2) .* Δ¹d))
end
function _quadratic_coeff(
    size_class_densities::Array{Float64,3},
    linear_extensions::AbstractMatrix{Float64},
    Δ¹d::AbstractMatrix{Float64}
)::Vector{Float64}
    return dropdims(
        sum((3 .* (size_class_densities .* (linear_extensions .^ 2) .* Δ¹d)), dims=(1, 2)),
        dims=(1, 2)
    )
end

function _linear_coeff(
    size_class_densities::AbstractMatrix{Float64},
    linear_extensions::AbstractMatrix{Float64},
    Δ²d::AbstractMatrix{Float64}
)::Float64
    return sum(3 .* (size_class_densities .* linear_extensions .* Δ²d))
end
function _linear_coeff(
    size_class_densities::Array{Float64,3},
    linear_extensions::AbstractMatrix{Float64},
    Δ²d::AbstractMatrix{Float64}
)::Vector{Float64}
    return dropdims(
        sum((3 .* (size_class_densities .* linear_extensions .* Δ²d)), dims=(1, 2)),
        dims=(1, 2)
    )
end

function _constant_coeff(
    size_class_densities::AbstractMatrix{Float64},
    adjusted_projected_cover::Float64,
    Δ³d::AbstractMatrix{Float64},
)::Float64
    return sum((size_class_densities .* Δ³d)) - ((12 / π) * adjusted_projected_cover)
end
function _constant_coeff(
    size_class_densities::Array{Float64,3},
    adjusted_projected_cover::Vector{Float64},
    Δ³d::AbstractMatrix{Float64},
)::Vector{Float64}
    return dropdims(
        sum(((size_class_densities .* Δ³d)), dims=(1, 2)),
        dims=(1, 2)
    ) .- ((12 / π) .* adjusted_projected_cover)
end
