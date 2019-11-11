function imshow(m::MapImages.MapImage, args...; kwargs...)
    return plt.imshow(m, extent=MapImages.grid_extent(m), args...; kwargs...)
end

function imshow(ax, m::MapImages.MapImage, args...; kwargs...)
    return ax.imshow(m, extent=MapImages.grid_extent(m), args...; kwargs...)
end
