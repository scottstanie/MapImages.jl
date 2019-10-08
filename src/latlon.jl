# TODO: to update just the width!
# julia> Sario.DemRsc(demrsc, width=4)
# Sario.DemRsc
#   width: Int64 4
#
# seems great for slicing/cropping
# Also: use @pack! to set just one field
# @pack! demrsc = newwidth


# TODO: make a macro for functions to always get definted for img and demrsc?
#
#
function find_overlap_idxs(asc_img, desc_img, asc_demrsc::DemRsc, desc_demrsc::DemRsc)
    left, right, bottom, top = intersection_corners(asc_demrsc, desc_demrsc)
    println(left, right, bottom, top)

    # asc_patch = asc_velo[32.3:30.71, -104.1:-102.31]
    # desc_patch = desc_velo[32.3:30.71, -104.1:-102.31]

    row1, col1 = nearest_pixel(asc_demrsc, top, left)
    row2, col2 = nearest_pixel(asc_demrsc, bottom, right)
    asc_idxs = (row1:row2, col1:col2)

    # asc_patch = asc_img[row1:row2, col1:col2]

    row1, col1 = nearest_pixel(desc_demrsc, top, left)
    row2, col2 = nearest_pixel(desc_demrsc, bottom, right)
    desc_idxs = (row1:row2, col1:col2)
    # desc_patch = desc_img[row1:row2, col1:col2]

    return asc_idxs, desc_idxs
    # return asc_patch, desc_patch
end


find_overlap_idxs(asc::MapImage, desc::MapImage) = find_overlap_idxs(asc.image, desc.image, asc.demrsc, desc.demrsc)

function find_overlap_patches(asc_fname::AbstractString, desc_fname::AbstractString=asc_fname, dset="velos/1"; mask=true)

    asc_img = MapImage(asc_fname, dset)
    desc_img = Sario.load(desc_fname, dset)

    asc_idxs, desc_idxs = find_overlap_idxs(asc_img, desc_img)

    a, d = asc_img[asc_idxs...], desc_img[desc_idxs...]
    if mask
        a, d = _mask_asc_desc(a, d)
    end
    return a, d
end

function _mask_asc_desc(a, d)
    m1 = a .== 0;
    m2 = d .== 0;
    mask = m1 .| m2
    a[mask] .= 0;
    d[mask] .= 0;
    return a, d
end

function find_overlap(asc::MapImage, desc::MapImage) 
    left, right, bottom, top = intersection_corners(asc, desc)

    row1, col1 = nearest_pixel(asc, top, left)
    row2, col2 = nearest_pixel(asc, bottom, right)
    asc_idxs = (row1:row2, col1:col2)
    asc_patch = asc_img[row1:row2, col1:col2]

    row1, col1 = nearest_pixel(desc, top, left)
    row2, col2 = nearest_pixel(desc, bottom, right)
    desc_idxs = (row1:row2, col1:col2)
    desc_patch = desc_img[row1:row2, col1:col2]

    return asc_patch, desc_patch
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
function nearest_pixel(demrsc::DemRsc, lat::AbstractFloat, lon::AbstractFloat)
    @unpack x_first, x_step, y_first, y_step = demrsc

    col_idx = 1 + (lon - x_first) / x_step
    # out_row_col[2] = _check_bounds(col_idx_arr, ncols)
    row_idx = 1 + (lat - y_first) / y_step
    # out_row_col[1] = _check_bounds(row_idx_arr, nrows)

    return Int.(round.((row_idx, col_idx)))
end

nearest_pixel(demrsc::DemRsc,
              lats::AbstractArray{<:AbstractFloat}, 
              lons::AbstractArray{<:AbstractFloat}) = [nearest_pixel(demrsc, lat, lon)
                                                       for (lat, lon) in zip(lats, lons)]


nearest_pixel(img::MapImage, lats, lons) = nearest_pixel(img.demrsc, lats, lons)
    

"""
Returns:
    tuple[float]: the boundaries of the intersection box of the 2 areas in order:
    (lon_left,lon_right,lat_bottom,lat_top)
"""
function intersection_corners(dem1::DemRsc, dem2::DemRsc)
    corners1 = grid_corners(dem1)
    corners2 = grid_corners(dem2)
    lons1, lats1 = zip(corners1...)
    lons2, lats2 = zip(corners2...)
    left = _max_min(lons1, lons2)
    right = _least_common(lons1, lons2)
    bottom = _max_min(lats1, lats2)
    top = _least_common(lats1, lats2)
    return left, right, bottom, top
end

intersection_corners(img1::MapImage, img2::MapImage) = intersection_corners(img1.demrsc, img2.demrsc)


_max_min(a, b) = max(minimum(a), minimum(b))
_least_common(a, b) = min(maximum(a), maximum(b))


function grid_corners(demrsc::DemRsc)
    """Takes sizes and spacing from .rsc info, finds corner points in (x, y) form

    Returns:
        list[tuple[float]]: the corners of the latlon grid in order:
        (top right, top left, bottom left, bottom right)
    """
    left, right, bot, top = grid_extent(demrsc)
    return [(right, top), (left, top), (left, bot), (right, bot)]
end
grid_corners(img::MapImage) = grid_corners(img.demrsc)


"""
Get the boundaries of a grid

Returns:
    tuple[float]: the boundaries of the latlon grid in order:
    (lon_left,lon_right,lat_bottom,lat_top)
"""
function grid_extent(demrsc::DemRsc)
    @unpack rows, cols, y_step, x_step, y_first, x_first = demrsc
    return (x_first, x_first .+ x_step * (cols - 1), y_first + y_step * (rows - 1), y_first)
end
grid_extent(img::MapImage) = grid_extent(img.demrsc)

"""Returns meshgrid-like arrays X, Y

julia> XX, YY = MapImages.grid(demrsc);

julia> XX[1:3, 1:3]
3×3 Array{Float64,2}:
 -104.0  -103.998  -103.997
 -104.0  -103.998  -103.997
 -104.0  -103.998  -103.997

julia> YY[1:3, 1:3]
3×3 Array{Float64,2}:
 31.9     31.9     31.9
 31.8983  31.8983  31.8983
 31.8967  31.8967  31.8967

julia> XX, YY = MapImages.grid(demrsc, sparse=true);
"""
function grid(demrsc::DemRsc; sparse=false)
    @unpack rows, cols, y_step, x_step, y_first, x_first = demrsc
    x = range(x_first, step=x_step, length=cols)
    y = range(y_first, step=y_step, length=rows)
    sparse && return (x, y)

    # return collect(reshape(y, :, 1) .* reshape(x, 1, :))
    return collect(ones(length(y), 1) .* reshape(x, 1, :)), collect(reshape(y, :, 1) .* ones(1, length(x)))
end
grid(img::MapImage; sparse=false) = grid(img.demrsc, sparse=sparse)


""" Takes the row, col of a pixel and finds its lat/lon """
function rowcol_to_latlon(row, col, demrsc::DemRsc)
    @unpack y_step, x_step, y_first, x_first = demrsc
    lat = y_first + (row - 1) * y_step
    lon = x_first + (col - 1) * x_step
    return lat, lon
end



"""Find the distance between two lat/lon points on Earth

Uses the haversine formula: https://en.wikipedia.org/wiki/Haversine_formula
so it does not account for the ellopsoidal Earth shape. Will be with about
0.5-1% of the correct value.

Notes: lats and lons are in degrees, and the values used for R Earth
(6373 km) are optimized for locations around 39 degrees from the equator

Reference: https://andrew.hedges.name/experiments/haversine/
"""
function latlon_to_dist(lat_lon_start, lat_lon_end, R=6378)
    lat1, lon1 = lat_lon_start
    lat2, lon2 = lat_lon_end
    dlat = deg2rad(lat2 - lat1)
    dlon = deg2rad(lon2 - lon1)
    lat1 = deg2rad(lat1)
    lat2 = deg2rad(lat2)
    a = (sin(dlat / 2)^2) + (cos(lat1) * cos(lat2) * sin(dlon / 2)^2)
    c = 2 * atan(sqrt(a), sqrt(1 - a))
    return R * c

end
