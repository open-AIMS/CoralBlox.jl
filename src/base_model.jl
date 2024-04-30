module base_model
using YAXArrays

"""
Total available area [m²]
"""
k_max()::Int64 = 10000

"""
Settlers base fraction for each species[m²]
"""
settlers_base()::Vector{Float64} = k_max() .* [0.1, 0.15]

"""
Diameter size class bins [m]
"""
bins()::Matrix{Float64} = [0.0 0.01 0.2 0.4 0.5; 0.0 0.01 0.2 0.4 0.5]

"""
Diameter Linear Extension [m/year]
"""
linear_extension()::Matrix{Float64} = [0.005 0.02 0.02 0.0; 0.004 0.017 0.023 0.0]

"""
Background Mortality rate [1/year]
"""
base_mortality_rate()::Matrix{Float64} = [0.02 0.01 0.01 0.005; 0.019 0.012 0.013 0.007]

"""
Fraction of k_max filled with each size class corals
"""
function init_cover()
    init_cover_fracs = [0.1 0.05 0.05 0.01; 0.15 0.07 0.03 0.02]
    n_species, n_bins = size(init_cover_fracs)
    axlist = (Dim{:species}(1:n_species), Dim{:size}(1:n_bins))
    rel_init_cover = YAXArray(axlist, init_cover_fracs)
    @assert (sum(rel_init_cover) <= 1) "Sum of all relative cover should be smaller than or equal to 1"

    # Initial coral cover
    return k_max() .* rel_init_cover
end

"""
Fraction of k_max available.
"""
function _k(cover::Array{Float64})::Float64
    sum(cover) > k_max() ? error("k-area should not be greater than k_max") : nothing
    return 1 - (sum(cover) / k_max())
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
function _diameter_factor(bins, s::Real, species::Real, di::Real, df::Real)
    res = (df^3 - di^3) / (bins[species, s+1]^3 - bins[species, s]^3)
    return res
end

"""
- `cover` : timesteps X species X size class
- `t` : Time step of that iteration
"""
function timestep_iteration(
    cover::AbstractArray{Float64},
    t::Int64,
    bins::Matrix{Float64},
    linear_extension::Matrix{Float64},
    base_mortality_rate::Matrix{Float64}
)
    n_species, n_bins = size(cover)[2:3]
    cover_t::Matrix{Float64} = zeros(Float64, n_species, n_bins)

    # Available area
    k::Float64 = _k(cover[t-1, :, :])

    # Settlers [area] for each species
    settlers::Vector{Float64} = settlers_base() .* log(1 + k)

    survival_rate::Matrix{Float64} = 1 .- base_mortality_rate
    growth_rate::Matrix{Float64} = linear_extension .* k

    small = 1
    medium = collect(2:(n_bins-1))
    large = n_bins

    species::Matrix{Int64} = [1 2; 1 2]

    m_growth_di::Matrix{Float64} = bins[:, medium] .+ growth_rate[:, medium]
    m_growth_df::Matrix{Float64} = bins[:, medium.+1]
    m_growth = _diameter_factor.([bins], [medium medium]', species, m_growth_di, m_growth_df)

    m_migrate_di = bins[:, medium]
    m_migrate_df = bins[:, medium] .+ growth_rate[:, medium.-1]
    m_migrate = _diameter_factor.([bins], [medium .- 1 medium .- 1]', species, m_migrate_di, m_migrate_df)

    l_migrate_factor = (3 .* growth_rate[:, large-1] .* bins[:, 4] .^ 2) ./ (bins[:, large] .^ 2 .- bins[:, large-1] .^ 2)

    cover_t[:, small] .= settlers
    @. cover_t[:, medium] = cover[t-1, :, medium] * survival_rate[:, medium] * m_growth +
                            cover[t-1, :, medium.-1] * survival_rate[:, medium.-1] * m_migrate
    @. cover_t[:, large] = cover[t-1, :, large] * survival_rate[:, large] +
                           cover[t-1, :, large-1] * survival_rate[:, large-1] .* l_migrate_factor

    @info "Cover after iteration $t-1: $(cover[t-1,:,:])"
    @info "K-area $k"
    return cover_t
end

function run_model()
    n_timesteps::Int64 = 100
    n_species::Int64 = 2
    n_bins::Int64 = 4
    cover::Array{Float64} = zeros(Float64, n_timesteps, n_species, n_bins)

    total_cover = zeros(n_timesteps, n_species)
    cover[1, :, :] = init_cover()
    total_cover[1, :] = sum(cover[1, :, :], dims=2)

    # axlist = (Dim{:species}(1:n_species), Dim{:size}(1:n_bins))
    _bins::Matrix{Float64} = bins()
    _linear_extension::Matrix{Float64} = linear_extension()
    _base_mortality_rate::Matrix{Float64} = base_mortality_rate()

    for t in 2:n_timesteps
        cover[t, :, :] .= timestep_iteration(cover, t, _bins, _linear_extension, _base_mortality_rate)
        total_cover[t, :] = sum(cover[t, :, :], dims=2)
    end

    return cover, total_cover
end
end
