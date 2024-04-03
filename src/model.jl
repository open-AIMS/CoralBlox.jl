# One location
# One species
# Three size classes

# Initial values
# Available area [m²]
k_max = 1000

# Settlers
r = k_max * 0.4

# Fraction of k_max filled with each size class corals []
rel_init_cover = [0.2, 0.2, 0.2, 0.01]

# Initial coral cover
init_cover = k_max .* rel_init_cover

# Diameter size class bins [m]
bins = [0.0, 0.01, 0.2, 0.4, 0.5]
n_bins = length(bins) - 1

# Diameter Linear Extension [m/year]
# Large corals don't expand
l = [0.005, 0.02, 0.02]

# Background Mortality rate [1/year]
bm = [0.02, 0.01, 0.01, 0.005]


"""
Fraction of k_max available.
"""
function _k(cover::Vector{Float64})
    return sum(cover) > k_max ? error("Problems") : 1 - (sum(cover) / k_max)
end

"""
- `s` : Size class
"""
function _bins_factor(s, di, df)
    return (df^3 - di^3) / (bins[s+1]^3 - bins[s]^3)
end

"""
- `cover` : timesteps X sizeclass Matrix
- `t` : Time step of that iteration
"""
function timestep_iteration(cover, t)
    cover_t::Vector{Float64} = zeros(Float64, n_bins)
    @info "Cover before iteration $t: $(cover[t-1, :])"
    # Available area
    k = _k(cover[t-1, :])

    # Settlers [area]
    ζ = r * log(1 + k)
    @info "k: $k"
    # Survival rate
    survivals = 1 .- bm
    #cover[t, :] .= cover[t, :] .* survivals

    # Increase in diameter
    D = l .* k

    small = 1
    medium = 2:(n_bins-1)
    large = n_bins

    cover_t[small] = ζ
    @. cover_t[medium] = cover[t-1, medium] * survivals[medium] * _bins_factor.(medium, bins[medium], (bins[medium.+1] - D[medium])) +
                         cover[t-1, medium.-1] * survivals[medium.-1] * _bins_factor.(medium .- 1, bins[medium] - D[medium.-1], bins[medium])
    cover_t[large] = cover[t-1, large] * survivals[large] +
                     cover[t-1, large-1] * survivals[large-1] * 3 * D[large-1] * bins[4]^2 / (bins[large]^2 - bins[large-1]^2)

    @info "Cover after iteration $t: $(cover_t)"
    return cover_t
end


# What to do with juveniles going to medium size class when there's no space to grow?

n_timesteps = 50
cover = zeros(n_timesteps, n_bins)
total_cover = zeros(n_timesteps)
cover[1, :] = init_cover
total_cover[1] = sum(cover[1, :])
for t in 2:n_timesteps
    cover[t, :] .= timestep_iteration(cover, t)
    total_cover[t] = sum(cover[t, :])
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
