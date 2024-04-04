using YAXArrays
using Plots
using Colors

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
function _diameter_factor(s::Real, species::Real, di::Real, df::Real)
    res = (df^3 - di^3) / (bins[species, s+1]^3 - bins[species, s]^3)
    return res
end

"""
- `cover` : timesteps X species X size class
- `t` : Time step of that iteration
"""
function timestep_iteration(cover, t)
    cover_t::Matrix{Float64} = zeros(Float64, n_species, n_bins)

    # Available area
    k = _k(cover[t-1, :, :])

    # Settlers [area]
    ζ = r * log(1 + k)

    # Survival rate
    survivals = 1 .- bm

    # Increase in diameter
    D = l .* k

    small = 1
    medium = collect(2:(n_bins-1))
    large = n_bins
    species = [1 2; 1 2]

    m_growth_di = bins[:, medium] .+ D[:, medium].data
    m_growth_df = bins[:, medium.+1]# .- D[:, medium].data
    m_growth = _diameter_factor.([medium medium]', species, m_growth_di, m_growth_df)

    m_migrate_di = bins[:, medium]
    m_migrate_df = bins[:, medium] .+ D[:, medium.-1].data
    m_migrate = _diameter_factor.([medium .- 1 medium .- 1]', species, m_migrate_di, m_migrate_df)

    l_migrate_factor = (3 .* D[:, large-1] .* bins[:, 4] .^ 2) ./ (bins[:, large] .^ 2 .- bins[:, large-1] .^ 2)

    cover_t[:, small] .= ζ
    @. cover_t[:, medium] = cover[t-1, :, medium] * survivals[:, medium] * m_growth +
                            cover[t-1, :, medium.-1] * survivals[:, medium.-1] * m_migrate
    @. cover_t[:, large] = cover[t-1, :, large] * survivals[:, large] +
                           cover[t-1, :, large-1] * survivals[:, large-1] .* l_migrate_factor

    @info "Cover after iteration $t-1: $(cover[t-1,:,:])"
    @info "K-area $k"
    return cover_t
end


# What to do with juveniles going to medium size class when there's no space to grow?

n_timesteps = 100
cover = zeros(n_timesteps, n_species, n_bins)
total_cover = zeros(n_timesteps, n_species)
cover[1, :, :] = init_cover
total_cover[1, :] = sum(cover[1, :, :], dims=2)
for t in 2:n_timesteps
    cover[t, :, :] .= timestep_iteration(cover, t)
    total_cover[t, :] = sum(cover[t, :, :], dims=2)
end

p_total = plot(
    collect(1:n_timesteps),
    [sum(total_cover[1:end, :], dims=2)],
    title="Total Coral Cover (k_area = $k_max)",
    label="Species 1",
    linewidth=3,
    xlabel="Timesteps",
    ylabel="Area",
    color=colorant"hsla(200, 100%, 50%, 1)",
    size=(500, 300)
)
savefig(p_total, "./figures/total_all_species_cover.png")

p_size_classes1 = plot(
    collect(1:n_timesteps),
    cover[1:end, 1, :],
    title="Species 1",
    label=["Small " "Medium 1" "Medium 2" "Large"],
    linewidth=2,
    xlabel="",
    ylabel="Area",
    color_palette=sequential_palette(200, 6)[end-4:end]
)

p_size_classes2 = plot(
    collect(1:n_timesteps),
    cover[:, 2, :],
    title="Species 2",
    label=["Small" "Medium 1" "Medium 2" "Large"],
    linewidth=2,
    xlabel="Timesteps",
    ylabel="Area",
    color_palette=sequential_palette(300, 6)[end-4:end]
)

p_each_species = plot(p_size_classes1, p_size_classes2, layout=(2, 1), size=(500, 600))
savefig(p_each_species, "./figures/each_species_cover.png")
