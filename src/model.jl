# One location
# One species
# Three size classes

# Initial values
# Available area [m²]
k_area = 1000

# Fraction of k_area filled with each size class corals []
rel_init_cover = [0.3, 0.2, 0.1]

# Diameter size class bins [m]
bins = [0.0, 0.01, 0.1, 0.5]

# Diameter Linear Extension [m/year]
Δd = [0.005, 0.02, 0.01]

# Background Mortality rate [1/year]
bm = [0.05, 0.03, 0.01]

init_cover = k_area .* rel_init_cover

function rel_available_space(cover::Vector{Float64})
    return max((1 - (sum(cover) / k_area)), 0.0)
end

# TODO check what happens in the limit available_space = 0
function timestep_iteration(cover_t)
    n_bins = length(cover_t)
    cover_t_1 = zeros(Float64, n_bins)

    # First apply mortality to all size classes
    @. cover_t_1 = cover_t - (cover_t * [0.05, 0.03, 0.01])
    @info "Cover t+1 with mortality: $cover_t_1"

    γ = rel_available_space(cover_t_1) .* Δd
    @info "Relative available space: $(rel_available_space(cover_t_1))"
    @info "Adjusted groth rate: $γ"

    # Cover loss to next size class
    # loss_cover = [1, 0.01, 0] .* cover_t_1
    loss_cover = zeros(Float64, n_bins)
    loss_cover[1] = cover_t_1[1]
    loss_cover[2] = ((bins[3] - (bins[3] - γ[2])^3) / (bins[3] - bins[2])) * cover_t_1[2]
    # @. loss_cover = (((bins[2:end] - γ)^3 - bins[2:end]^3) / (bins[2:end]^3 - bins[1:end-1]^3)) * cover_t_1
    loss_cover[3] = 0.0

    # Cover income from previous size class
    income_cover = zeros(Float64, n_bins)
    income_cover[1] = 0.0
    income_cover[2] = (cover_t_1[1] / (bins[2]^3 - bins[1]^3)) * ((bins[2] + γ[1])^3 - (bins[1] + γ[1])^3)
    income_cover[3] = (cover_t_1[2] / (bins[3]^3 - bins[2]^3)) * ((bins[3] + γ[2])^3 - (bins[3])^3)
    @info "Income factor 2: $(((bins[2] + γ[1])^3 - (bins[1] + γ[1])^3)/(bins[2]^3 - bins[1]^3))"

    # Cover internal growth
    internal_growth = zeros(Float64, n_bins)
    internal_growth[1] = 0
    internal_growth[2] = (cover_t_1[2] / (bins[3]^3 - bins[2]^3)) * ((bins[3])^3 - (bins[2] + γ[2])^3)
    internal_growth[3] = (cover_t_1[3] / (bins[4]^3 - bins[3]^3)) * ((bins[4] + γ[3])^3 - (bins[3] + γ[3])^3)

    # Settlers
    settlers = [200, 0, 0] .* rel_available_space(cover_t_1)

    @info "Income cover: $income_cover"
    @info "Internal growth: $internal_growth"
    # @info "Loss cover: $loss_cover"
    @info "Settlers: $settlers"

    @. cover_t_1 = income_cover + internal_growth + settlers

    return cover_t_1
end


# What to do with juveniles going to medium size class when there's no space to grow?

n_timesteps = 50
covers = zeros(3, n_timesteps)
covers[:, 1] = init_cover
for t in 2:n_timesteps
    covers[:, t] = timestep_iteration(covers[:, t-1])
end

using Plots

plot(collect(1:n_timesteps), [covers[1, :], covers[2, :], covers[3, :]], title="Teste", label=["small" "medium" "large"], linewidth=3)
