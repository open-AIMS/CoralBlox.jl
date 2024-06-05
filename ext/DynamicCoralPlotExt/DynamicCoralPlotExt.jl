module DynamicCoralPlotExt

using DynamicCoralCoverModel
import DynamicCoralCoverModel.blocks_model: CoverBlock
import DynamicCoralCoverModel.blocks_model: SizeClass
import DynamicCoralCoverModel.blocks_model: Δinterval
import DynamicCoralCoverModel.blocks_model: interval_lower_bound

using Makie, Makie.GeometryBasics

function rectangle(w, h, x, y)
    return x .+ [0, w, w, 0], y .+ [0, 0, h, h]
end

function rectangles!(sc, ind)
    for cv in sc.cover_blocks
        width = Δinterval(cv.interval)
        height = cv.diameter_density * (Δinterval(cv.interval))
        poly!(Rect(interval_lower_bound(cv), 0.0, width, log(10, height)), color=ind, strokewidth=0.0, colormap=:PuBuGn_6, colorrange=(1, 6))
    end
end

"""
This function will plot one figure for each timestep
"""
function DynamicCoralCoverModel.Plot.plot_size_class(size_classes::Matrix{SizeClass}, timestep::Int64)
    f = Figure(; size=(1600, 1600))
    n_species, n_bins = size(size_classes)
    x = [1, 1, 2, 2, 3, 3]
    y = [1, 2, 1, 2, 1, 2]
    for species in 1:n_species
        Axis(f[x[species], y[species]], title="spec: $(species), timestep = $timestep")

        for (size, sc) in enumerate(size_classes[species, 2:end])
            rectangles!(sc, size)
            # _x = vcat([[cv.interval[1], cv.interval[2]] for cv in sc.cover_blocks]...)
            # _y = vcat([fill(cv.diameter_density * (Δinterval(cv.interval)), 2) for cv in sc.cover_blocks]...)
            # n_corals += sum(_y) / 2
            # bar!(
            #     _x,
            #     _y,
            #     yscale=:log10,
            #     color=colors[size],
            #     fill=(0, 1, colors[size]),
            #     label="Size: $size",
            #     leg_title="SizeClass"
            # )
        end

        # plot!(
        #     [coral_spec.bins[1], coral_spec.bins[end]],
        #     [n_corals, n_corals],
        #     yscale=:log10,
        #     line=(2, :solid, :black)
        # )
    end

    save("./figures/gif/t$(timestep).png", f)
    return nothing
end

# function DynamicCoralCoverModel.Plot.plot_taxa(cover::Array{Float64, 3})
#     cover_t = dropdims(sum(cover, dims=3), dims=3) ./ 48671.0938
#     colors = [:red, :green, :blue, :purple, :yellow]
#
#     f = Figure()
#     Axis(f[1, 1], title="Relative Taxa Cover")
#     p1 = CairoMakie.lines!(1:100, cover_t[:, 1], color=colors[1])
#     p2 = CairoMakie.lines!(1:100, cover_t[:, 2], color=colors[2])
#     p3 = CairoMakie.lines!(1:100, cover_t[:, 3], color=colors[3])
#     p4 = CairoMakie.lines!(1:100, cover_t[:, 4], color=colors[4])
#     p5 = CairoMakie.lines!(1:100, cover_t[:, 5], color=colors[5])
#
#     Legend(f[1, 2],
#         [p1, p2, p3, p4, p5],
#         [
#            "Tabular Acropora",
#             "Corymbose Acropora",
#             "Corymbose non-Acropora",
#             "Small Massive",
#             "Large Massive"
#         ]
#     )
#
#     return f
# end
#
# function DynamicCoralCoverModel.Plot.plot_size_classes(cover, species, name)
#
#     cover = cover ./ 48671.0938
#
#     colors_alp = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#
#     f = Figure()
#     Axis(f[1, 1], title=name*" size classes")
#
#     p1 = CairoMakie.lines!(1:100, cover[:, species, 1], color=1, colormap=:PuBuGn_7, colorrange=(1, 7))
#     p2 = CairoMakie.lines!(1:100, cover[:, species, 2], color=2, colormap=:PuBuGn_7, colorrange=(1, 7))
#     p3 = CairoMakie.lines!(1:100, cover[:, species, 3], color=3, colormap=:PuBuGn_7, colorrange=(1, 7))
#     p4 = CairoMakie.lines!(1:100, cover[:, species, 4], color=4, colormap=:PuBuGn_7, colorrange=(1, 7))
#     p5 = CairoMakie.lines!(1:100, cover[:, species, 5], color=5, colormap=:PuBuGn_7, colorrange=(1, 7))
#     p6 = CairoMakie.lines!(1:100, cover[:, species, 6], color=6, colormap=:PuBuGn_7, colorrange=(1, 7))
#     p7 = CairoMakie.lines!(1:100, cover[:, species, 7], color=7, colormap=:PuBuGn_7, colorrange=(1, 7))
#
#     Legend(f[1, 2],
#         [p1, p2, p3, p4, p5, p6, p7],
#         [
#             "1",
#             "2",
#             "3",
#             "4",
#             "5",
#             "6",
#             "7"
#         ]
#     )
#     return f
# end
#
# function DynamicCoralCoverModel.Plot.plot_total_cover(cover)
#     total_cover = dropdims(sum(cover, dims=(2, 3)), dims=(2, 3)) ./ 48671.0938
#
#     f = Figure()
#     Axis(f[1, 1], title="Relative Cover")
#     CairoMakie.lines!(f[1, 1], 1:100, total_cover, color=:black)
#
#     return f
# end


end
