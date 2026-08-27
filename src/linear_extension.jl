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
    LinearExtensionCache(bin_edges)

Pre-computed Δd matrices derived from `bin_edges`. Construct once per domain/scenario
and pass to `linear_extension_scale_factors` (the spatial competition factor) to avoid
recomputing these constant matrices on every timestep.
"""
struct LinearExtensionCache
    Δ¹d::Matrix{Float64}
    Δ²d::Matrix{Float64}
    Δ³d::Matrix{Float64}
end
function LinearExtensionCache(bin_edges::AbstractMatrix{Float64})
    tgt = @view bin_edges[:, 1:(end - 1)]
    return LinearExtensionCache(Δd(tgt, 1), Δd(tgt, 2), Δd(tgt, 3))
end

# Private: computation with pre-computed Δd matrices and target_linear_extensions view.
# All public overloads delegate here after computing Δd once.
#
# This computes the spatial competition factor: the growth penalty applied uniformly to
# all linear extensions to account for corals competing for a shared, limited habitable
# area (see the "How does the linear extension scale factor work?" section of the README).
#
# density_i = (12/π) * C_i / Δ³_i is substituted into all coefficient expressions:
#   projected_cover = Σ C*(3*le²*Δ1/Δ3 + 3*le*Δ2/Δ3 + 1)
#   a = (36/π) * Σ C*le²*Δ1/Δ3
#   b = (36/π) * Σ C*le*Δ2/Δ3
#   c = (12/π) * (total_cover - adjusted_projected_cover)
# All three accumulators are computed in a single pass to avoid any intermediate allocations.
function _linear_extension_scale_factors_core(
    C_cover_t::AbstractMatrix{Float64},
    habitable_area::Float64,
    target_linear_extensions::AbstractMatrix{Float64},
    max_projected_cover::Float64,
    Δ¹d::AbstractMatrix{Float64},
    Δ²d::AbstractMatrix{Float64},
    Δ³d::AbstractMatrix{Float64},
)::Float64
    target_C_cover_t = @view C_cover_t[:, 1:(end - 1)]
    total_cover::Float64 = sum(target_C_cover_t)
    non_target_total_cover = sum(C_cover_t[:, end])

    projected_cover::Float64 = 0.0
    a_acc::Float64 = 0.0
    b_acc::Float64 = 0.0
    @inbounds for j ∈ axes(target_C_cover_t, 2), i ∈ axes(target_C_cover_t, 1)
        c_ij = target_C_cover_t[i, j]
        le = target_linear_extensions[i, j]
        Δ1 = Δ¹d[i, j]
        Δ2 = Δ²d[i, j]
        Δ3 = Δ³d[i, j]
        c_over_Δ3 = c_ij / Δ3
        projected_cover += c_ij + c_over_Δ3 * (3 * le^2 * Δ1 + 3 * le * Δ2)
        a_acc += c_over_Δ3 * le^2 * Δ1
        b_acc += c_over_Δ3 * le * Δ2
    end

    a::Float64 = (36 / π) * a_acc
    b::Float64 = (36 / π) * b_acc

    adjusted_projected_cover::Float64 = _adjusted_projected_cover(
        total_cover,
        projected_cover,
        max_projected_cover - non_target_total_cover,
        habitable_area - non_target_total_cover,
    )

    c::Float64 = (12 / π) * (total_cover - adjusted_projected_cover)

    return (sqrt((b^2) - (4 * a * c)) - b) / (2 * a)
end

# Public scalar: computes Δd once then delegates to core.
"""
    linear_extension_scale_factors(C_cover_t::AbstractArray{Float64,3}, loc_habitable_areas::AbstractVector{Float64}, linear_extensions::AbstractMatrix{Float64}, bin_edges::AbstractMatrix{Float64}, max_projected_cover::AbstractVector{Float64})

Computes the spatial competition factor: a growth penalty, uniform across `FunctionalGroup`s
and `SizeClass`es, that captures how growth slows down purely because corals are competing
for a shared, limited habitable area. Multiplying `linear_extensions` by this factor gives
the adjusted linear extension. It assumes the last size class of each functional group doesn't
grow, so those size classes are excluded from this calculation to prevent a
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
    target_bin_edges = @view bin_edges[:, 1:(end - 1)]
    Δ¹d = Δd(target_bin_edges, 1)
    Δ²d = Δd(target_bin_edges, 2)
    Δ³d = Δd(target_bin_edges, 3)
    return _linear_extension_scale_factors_core(
        C_cover_t,
        habitable_area,
        @view(linear_extensions[:, 1:(end - 1)]),
        max_projected_cover,
        Δ¹d,
        Δ²d,
        Δ³d,
    )
end

# Public vector: computes Δd once for the whole batch, then calls core per location.
function linear_extension_scale_factors(
    C_cover_t::AbstractArray{Float64, 3},
    loc_habitable_areas::AbstractVector{Float64},
    linear_extensions::AbstractMatrix{Float64},
    bin_edges::AbstractMatrix{Float64},
    max_projected_cover::AbstractVector{Float64},
)::AbstractVector{Float64}
    return linear_extension_scale_factors(
        C_cover_t,
        loc_habitable_areas,
        linear_extensions,
        LinearExtensionCache(bin_edges),
        max_projected_cover,
    )
end

# Public vector with a LinearExtensionCache — for callers that pre-compute the cache
# once before a timestep loop to avoid recomputing constant Δd matrices each call.
function linear_extension_scale_factors(
    C_cover_t::AbstractArray{Float64, 3},
    loc_habitable_areas::AbstractVector{Float64},
    linear_extensions::AbstractMatrix{Float64},
    cache::LinearExtensionCache,
    max_projected_cover::AbstractVector{Float64},
)::AbstractVector{Float64}
    n = size(C_cover_t, 3)
    result = Vector{Float64}(undef, n)
    target_lin_ext = @view linear_extensions[:, 1:(end - 1)]
    @views for i ∈ 1:n
        result[i] = _linear_extension_scale_factors_core(
            C_cover_t[:, :, i],
            loc_habitable_areas[i],
            target_lin_ext,
            max_projected_cover[i],
            cache.Δ¹d,
            cache.Δ²d,
            cache.Δ³d,
        )
    end
    return result
end

@inline @views Δd(bin_edges::AbstractMatrix{Float64}, n::Int64)::Matrix{Float64} =
    @. ((bin_edges[:, 2:end]^n) - (bin_edges[:, 1:(end - 1)]^n))

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
