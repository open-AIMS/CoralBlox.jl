"""
    max_projected_cover(linear_extensions::Matrix{Float64}, bin_edges::Matrix{Float64}, habitable_area::Float64)::Float64
    max_projected_cover(linear_extensions::Matrix{Float64}, bin_edges::Matrix{Float64}, habitable_areas::Vector{Float64})

Finds the maximum projected cover for each location considering that all habitable area is
concentrated on a single size class and functional group and that it grows by that
functional group's linear extension. To choose which combination of functional group and
size class all combinations are tested.

# Arguments
- `linear_extensions` : Functional Groups x Size Classes
- `bin_edges` : Functional Groups x Bin Edges
- `habitable_areas` : Habitable area
"""
function max_projected_cover(
    linear_extensions::Matrix{Float64}, bin_edges::Matrix{Float64}, habitable_area::Float64
)::Float64
    return _diameter_coef(linear_extensions, bin_edges) * habitable_area
end
function max_projected_cover(
    linear_extensions::Matrix{Float64},
    bin_edges::Matrix{Float64},
    habitable_areas::Vector{Float64},
)::Vector{Float64}
    return [_diameter_coef(linear_extensions, bin_edges)] .* habitable_areas
end

function _diameter_coef(
    linear_extensions::Matrix{Float64}, bin_edges::Matrix{Float64}
)::Float64
    size_class_coefs =
        ((bin_edges[:, 2:(end - 1)] .+ linear_extensions[:, 1:(end - 1)]) .^ 2) ./
        (bin_edges[:, 2:(end - 1)] .^ 2)
    return findmax(size_class_coefs)[1]
end

"""
    adjusted_linear_extension(C_cover_t::AbstractArray{Float64,3}, loc_habitable_areas::AbstractVector{Float64}, linear_extensions::AbstractMatrix{Float64}, bin_edges::AbstractMatrix{Float64}, max_projected_cover::AbstractVector{Float64})

Adjusted linear extension. It assumes the last functional group doesn't grow. Therefore,
the last size class of each functional group are excluded from this calculation to prevent a
rounding error that sometimes occurs when almost all of the cover is concentrated in the
last size class of some functional groups and the growth (correspondent to the remaining
size classes) is marginal and, because of rounding, sometimes makes the projected_cover be
slightly larger than loc_C_cover for some locs.

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
    max_projected_cover::Float64,
)::Float64
    # Target here refers to all size classes except the last one of each functional group
    # since these don't grow
    target_linear_extensions = @view linear_extensions[:, 1:(end - 1)]
    target_bin_edges = @view bin_edges[:, 1:(end - 1)]
    target_C_cover_t = @view C_cover_t[:, 1:(end - 1)]

    total_cover::Float64 = sum(target_C_cover_t)
    non_target_total_cover = sum(C_cover_t[:, end])

    # Average density for each functional group and size class
    size_class_densities::Matrix{Float64} = _size_class_densities(
        target_C_cover_t, target_bin_edges
    )

    # Projected coral cover at t+1
    projected_cover::Float64 = _projected_cover(
        size_class_densities, target_linear_extensions, target_bin_edges
    )

    adjusted_projected_cover::Float64 = _adjusted_projected_cover(
        total_cover,
        projected_cover,
        max_projected_cover - non_target_total_cover,
        habitable_area - non_target_total_cover,
    )

    # Solve quadratic equation
    a::Float64 = _quadratic_coeff(
        size_class_densities, target_linear_extensions, Δd(target_bin_edges, 1)
    )
    b::Float64 = _linear_coeff(
        size_class_densities, target_linear_extensions, Δd(target_bin_edges, 2)
    )
    c::Float64 = _constant_coeff(
        size_class_densities, adjusted_projected_cover, Δd(target_bin_edges, 3)
    )

    return (sqrt((b^2) - (4 * a * c)) - (b)) / (2 * a)
end
function linear_extension_scale_factors(
    C_cover_t::AbstractArray{Float64, 3},
    loc_habitable_areas::AbstractVector{Float64},
    linear_extensions::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64},
    max_projected_cover::AbstractVector{Float64},
)::AbstractVector{Float64}
    n = size(C_cover_t, 3)
    result = Vector{Float64}(undef, n)
    @views for i ∈ 1:n
        result[i] = linear_extension_scale_factors(
            C_cover_t[:, :, i],
            loc_habitable_areas[i],
            linear_extensions,
            bin_edges,
            max_projected_cover[i],
        )
    end
    return result
end

Δd(bin_edges::AbstractMatrix{Float64}, n::Int64)::Matrix{Float64} =
    ((bin_edges[:, 2:end] .^ n) .- (bin_edges[:, 1:(end - 1)] .^ n))

_size_class_densities(
    C_cover_t::AbstractMatrix{Float64}, bin_edges::AbstractMatrix{Float64}
)::Matrix{Float64} = (12 / π) .* (C_cover_t ./ Δd(bin_edges, 3))
_size_class_densities(
    C_cover_t::AbstractArray{Float64, 3}, bin_edges::AbstractMatrix{Float64}
)::Array{Float64, 3} = (12 / π) .* (C_cover_t ./ Δd(bin_edges, 3))

function _projected_cover(
    size_class_densities::AbstractMatrix{Float64},
    linear_extensions::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64},
)::Float64
    return sum(
        (π / 12) .* size_class_densities .* (
            (3 .* (linear_extensions .^ 2) .* Δd(bin_edges, 1)) .+
            (3 .* linear_extensions .* Δd(bin_edges, 2)) .+ Δd(bin_edges, 3)
        ),
    )
end
function _projected_cover(
    size_class_densities::AbstractArray{Float64, 3},
    linear_extensions::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64},
)::Vector{Float64}
    return dropdims(
        sum(
            (
                (π / 12) .* size_class_densities .* (
                    (3 .* (linear_extensions .^ 2) .* Δd(bin_edges, 1)) .+
                    (3 .* linear_extensions .* Δd(bin_edges, 2)) .+ Δd(bin_edges, 3)
                )
            );
            dims=(2, 1),
        );
        dims=(2, 1),
    )
end

function _adjusted_projected_cover(
    total_cover::Float64,
    projected_cover::Float64,
    max_projected_cover::Float64,
    habitable_area::Float64,
)::Float64
    return total_cover + (
        (projected_cover - total_cover) *
        ((habitable_area - total_cover) / (max_projected_cover - total_cover))
    )
end
function _adjusted_projected_cover(
    loc_cover::AbstractVector{Float64},
    projected_cover::AbstractVector{Float64},
    max_projected_cover::AbstractVector{Float64},
    loc_habitable_areas::AbstractVector{Float64},
)::Vector{Float64}
    return loc_cover .+ (
        (projected_cover .- loc_cover) .*
        ((loc_habitable_areas .- loc_cover) ./ (max_projected_cover .- loc_cover))
    )
end

function _quadratic_coeff(
    size_class_densities::AbstractMatrix{Float64},
    linear_extensions::AbstractMatrix{Float64},
    Δ¹d::AbstractMatrix{Float64},
)::Float64
    return sum(3 .* (size_class_densities .* (linear_extensions .^ 2) .* Δ¹d))
end
function _quadratic_coeff(
    size_class_densities::Array{Float64, 3},
    linear_extensions::AbstractMatrix{Float64},
    Δ¹d::AbstractMatrix{Float64},
)::Vector{Float64}
    return dropdims(
        sum((3 .* (size_class_densities .* (linear_extensions .^ 2) .* Δ¹d)); dims=(1, 2));
        dims=(1, 2),
    )
end

function _linear_coeff(
    size_class_densities::AbstractMatrix{Float64},
    linear_extensions::AbstractMatrix{Float64},
    Δ²d::AbstractMatrix{Float64},
)::Float64
    return sum(3 .* (size_class_densities .* linear_extensions .* Δ²d))
end
function _linear_coeff(
    size_class_densities::Array{Float64, 3},
    linear_extensions::AbstractMatrix{Float64},
    Δ²d::AbstractMatrix{Float64},
)::Vector{Float64}
    return dropdims(
        sum((3 .* (size_class_densities .* linear_extensions .* Δ²d)); dims=(1, 2));
        dims=(1, 2),
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
    size_class_densities::Array{Float64, 3},
    adjusted_projected_cover::Vector{Float64},
    Δ³d::AbstractMatrix{Float64},
)::Vector{Float64}
    return dropdims(sum(((size_class_densities .* Δ³d)); dims=(1, 2)); dims=(1, 2)) .-
           ((12 / π) .* adjusted_projected_cover)
end
