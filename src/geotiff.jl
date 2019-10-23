function create_geotiff(img_filename, outfile=replace(img_filename, ".png" => ".tif");
                        demrsc=nothing, demrscfile=nothing)

    if isnothing(demrsc)
        demrsc = isnothing(demrscfile) ? Sario.load(Sario.find_rsc_file(".")) : Sario.load(demrscfile)
    end
    # gdalwarp -s_srs EPSG:4326 -t_srs EPSG:3857 out1.tif out2.tif
    # ullr means upper left, lower right (lon, lat) points
    north, south, east, west = nsew(demrsc)
    cmd = `gdal_translate -of GTiff -a_nodata 0 -a_srs EPSG:4326 -a_ullr $west $north $east $south $img_filename $outfile`
    println("Running:")
    println(cmd)
    run(cmd)
end

"""return tuple of (north, south, east, west) from demrsc data"""
function nsew(demrsc)
    @unpack x_first, y_first, x_step, y_step, width, file_length = demrsc
    north = y_first
    west = x_first
    east = west + width * x_step
    south = north + file_length * y_step
    return north, south, east, west
end
