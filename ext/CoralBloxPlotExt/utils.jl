function grid_positions(n_ax::Int)::Vector{Vector{Int}}
    n_cols = Int(floor(sqrt(n_ax)))
    n_rows = ceil(Int, n_ax / n_cols)
    return [
        [i, j] for i in 1:n_rows, j in 1:n_cols
    ][1:n_ax]
end
