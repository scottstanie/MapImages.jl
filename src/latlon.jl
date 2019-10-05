

function find_overlap_idxs(asc_img, desc_img, asc_dem_rsc, desc_dem_rsc)
    left, right, bottom, top = intersection_corners(asc_dem_rsc, desc_dem_rsc)
    println(left, right, bottom, top)

    # asc_patch = asc_velo[32.3:30.71, -104.1:-102.31]
    # desc_patch = desc_velo[32.3:30.71, -104.1:-102.31]

    row1, col1 = nearest_pixel(asc_dem_rsc, top, left)
    row2, col2 = nearest_pixel(asc_dem_rsc, bottom, right)
    asc_idxs = (row1:row2, col1:col2)

    # asc_patch = asc_img[row1:row2, col1:col2]

    row1, col1 = nearest_pixel(desc_dem_rsc, top, left)
    row2, col2 = nearest_pixel(desc_dem_rsc, bottom, right)
    desc_idxs = (row1:row2, col1:col2)
    # desc_patch = desc_img[row1:row2, col1:col2]

    return asc_idxs, desc_idxs
    # return asc_patch, desc_patch
end

function find_overlap_patches(asc_fname::AbstractString, desc_fname::AbstractString=asc_fname, dset="velos/1")

    asc_img = Sario.load(asc_fname, dset_name=dset)
    desc_img = Sario.load(desc_fname, dset_name=dset)

    asc_dem_rsc = Sario.load_dem_from_h5(asc_fname)
    desc_dem_rsc = Sario.load_dem_from_h5(desc_fname)
    asc_idxs, desc_idxs = find_overlap_idxs(asc_img, desc_img, asc_dem_rsc, desc_dem_rsc)

    a, d = asc_img[asc_idxs...], desc_img[desc_idxs...]
    m1 = a .== 0;
    m2 = d .== 0;
    mask = m1 .| m2
    a[mask] .= 0;
    d[mask] .= 0;
    return a, d
end

# function _check_bounds(idx_arr, bound)
#     int_idxs = Int.(round.(idx_arr))
#     bad_idxs = int_idxs .< 0 .| int_idxs .>= bound
# 
#     if any(bad_idxs)
#         # Need to check for single numbers, shape ()
#         if int_idxs.shape:
#             # Replaces locations of bad_idxs with none
#             int_idxs = findall(bad_idxs, None, int_idxs)
#         else:
#             int_idxs = None
#         end
#     end
# 
#     return int_idxs
# end

"""Find the nearest row, col to a given lat and/or lon"""
function nearest_pixel(dem_rsc, lat::AbstractFloat, lon::AbstractFloat)
    col_idx = 1 + (lon - dem_rsc["x_first"]) / dem_rsc["x_step"]
    # out_row_col[2] = _check_bounds(col_idx_arr, ncols)
    row_idx = 1 + (lat - dem_rsc["y_first"]) / dem_rsc["y_step"]
    # out_row_col[1] = _check_bounds(row_idx_arr, nrows)

    return Int.(round.((row_idx, col_idx)))
end
nearest_pixel(dem_rsc, lats::AbstractArray{AbstractFloat}, 
              lons::AbstractArray{AbstractFloat}) = [nearest_pixel(dem_rsc, lat, lon)
                                                     for (lat, lon) in zip(lats, lons)]

_max_min(a, b) = max(minimum(a), minimum(b))
_least_common(a, b) = min(maximum(a), maximum(b))

# TODO: switch all dems to symbols??
_symdict(d) = Dict(Symbol(k) => v for (k, v) in d)
function intersection_corners(dem1, dem2)
    """
    Returns:
        tuple[float]: the boundaries of the intersection box of the 2 areas in order:
        (lon_left,lon_right,lat_bottom,lat_top)
    """
    dem1 = _symdict(dem1)
    dem2 = _symdict(dem2)
    corners1 = grid_corners(;dem1...)
    corners2 = grid_corners(;dem2...)
    lons1, lats1 = zip(corners1...)
    lons2, lats2 = zip(corners2...)
    left = _max_min(lons1, lons2)
    right = _least_common(lons1, lons2)
    bottom = _max_min(lats1, lats2)
    top = _least_common(lats1, lats2)
    return left, right, bottom, top
end


function grid_corners(; kwargs...)
    """Takes sizes and spacing from .rsc info, finds corner points in (x, y) form

    Returns:
        list[tuple[float]]: the corners of the latlon grid in order:
        (top right, top left, bottom left, bottom right)
    """
    left, right, bot, top = grid_extent(;kwargs...)
    return [(right, top), (left, top), (left, bot), (right, bot)]
end


"""
Returns:
    tuple[float]: the boundaries of the latlon grid in order:
    (lon_left,lon_right,lat_bottom,lat_top)
"""
function grid_extent(;rows=nothing,
                     cols=nothing,
                     y_step=nothing,
                     x_step=nothing,
                     y_first=nothing,
                     x_first=nothing,
                     file_length=nothing,
                     width=nothing,
                     kwargs...)
    rows = isnothing(rows) ? file_length : rows
    cols = isnothing(cols) ? width : cols
    @show cols
    return (x_first, x_first .+ x_step * (cols - 1), y_first + y_step * (rows - 1), y_first)
end

