using YAXArrays

n_species = 2
n_bins = 4
axlist = (Dim{:species}(1:n_species), Dim{:size}(1:n_bins))

# Available area [m²]
k_max = 10000

# Settlers
r = k_max .* [0.1, 0.15]

# Fraction of k_max filled with each size class corals []
rel_init_cover = YAXArray(axlist, [0.1 0.05 0.05 0.01; 0.15 0.07 0.03 0.02])
@assert (sum(rel_init_cover) <= 1) "Sum of all relative cover should be smaller than or equal to 1"

# Initial coral cover
init_cover = k_max .* rel_init_cover

# Diameter size class bins [m]
bins = [0.0 0.01 0.2 0.4 0.5; 0.0 0.01 0.2 0.4 0.5]

# Diameter Linear Extension [m/year]
l = YAXArray(axlist, [0.005 0.02 0.02 0.0; 0.004 0.017 0.023 0.0])

# Background Mortality rate [1/year]
bm = YAXArray(axlist, [0.02 0.01 0.01 0.005; 0.019 0.012 0.013 0.007])

"""
Fraction of k_max available.
"""
function _k(cover::Array{Float64})
    return sum(cover) > k_max ? error("k-area should not be greater than k_max") : 1 - (sum(cover) / k_max)
end

"""
Diameter proportionality factor of size class `s` corals with diameter between `di` and
`df`. This factor, multiplied by the total cover `C_s` of size class `s` corals gives the
area of the corals within that size class with diameter between `di` and `df`.

# Arguments
- `s` : Size class
- `di`: Initial diameter
- `df`: Final diameter
"""
function _diameter_factor(s, di, df)
    @info "s = $s ; di = $di ; df = $df"
    @info (df^3 - di^3) ./ (bins[:, s+1] .^ 3 .- bins[:, s] .^ 3)
    return (df^3 - di^3) ./ (bins[:, s+1] .^ 3 .- bins[:, s] .^ 3)
end

"""
- `cover` : timesteps X species X size class
- `t` : Time step of that iteration
"""
function timestep_iteration(cover, t)
    cover_t::Matrix{Float64} = zeros(Float64, n_species, n_bins)

    # Available area
    k = _k(cover[t-1, :, :])
    k = _k(cover[1, :, :])

    # Settlers [area]
    ζ = r * log(1 + k)

    # Survival rate
    survivals = 1 .- bm

    # Increase in diameter
    D = l .* k

    small = 1
    medium = collect(2:(n_bins-1))
    large = n_bins

    cover_t[:, small] .= ζ
    @. cover_t[:, medium] = cover[t-1, :, medium] * survivals[:, medium] * _diameter_factor.([medium medium], bins[:, medium]', (bins[:, medium.+1] - D[:, medium])') +
                            cover[t-1, :, medium.-1] * survivals[:, medium.-1] * _diameter_factor.(medium .- 1, bins[:, medium] - D[:, medium.-1], bins[:, medium])
    @. cover_t[:, large] = cover[t-1, :, large] * survivals[:, large] +
                           cover[t-1, :, large-1] * survivals[:, large-1] .* (3 .* D[:, large-1] .* bins[:, 4]^2) ./ (bins[:, large]^2 .- bins[:, large-1]^2)

    @info "Cover after iteration $t: $(cover_t)"
    return cover_t
end


# What to do with juveniles going to medium size class when there's no space to grow?

n_timesteps = 50
cover = zeros(n_timesteps, n_species, n_bins)
total_cover = zeros(n_timesteps, n_species)
cover[1, :, :] = init_cover
total_cover[1, :] = sum(cover[1, :, :], dims=2)
for t in 2:n_timesteps
    cover[t, :, :] .= timestep_iteration(cover, t)
    total_cover[t, :] = sum(cover[t, :, :], dims=2)
end

using Plots
using Colors

p_total = plot(
    collect(1:n_timesteps),
    [total_cover],
    title="Total Coral Cover (k_area = $k_max)",
    label="Species 1",
    linewidth=3,
    xlabel="Timesteps",
    ylabel="Area",
    color=colorant"hsla(200, 100%, 50%, 1)"
)
savefig(p_total, "./figures/total_cover.png")

p_size_classes = plot(
    collect(1:n_timesteps),
    cover,
    title="Coral Cover by size class (k_area = $k_max)",
    label=["Small" "Medium 1" "Medium 2" "Large"],
    linewidth=2,
    xlabel="Timesteps",
    ylabel="Area",
    color_palette=sequential_palette(200, 6)[end-4:end]
)
savefig(p_size_classes, "./figures/size_classes.png")
