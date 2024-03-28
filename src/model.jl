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
Δd = [0.01, 0.03, 0.02]

# Background Mortality rate [1/year]
bm = [0.05, 0.03, 0.01]

init_cover = k_area .* rel_init_cover

function timestep_iteration(init_cover)
    # First apply mortality to all size classes
    init_cover .-= init_cover .* [0.05, 0.03, 0.01]


    # Cover income from previous size class
    # TODO
    income_cover = [0, 0.01, 0.01] .* init_cover

    # Cover loss to next size class
    # TODO
    loss_cover = [1, 0.01, 0] .* init_cover

    # Cover internal growth
    # TODO
    internal_growth = [0, 0.02, 0.02] .* init_cover

    # Settlers
    settlers = [0.1 * init_cover[3], 0, 0]

    @. init_cover += (income_cover + internal_growth - loss_cover + settlers)

    return init_cover
end
