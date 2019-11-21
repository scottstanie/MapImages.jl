var documenterSearchIndex = {"docs":
[{"location":"#MapImages.jl-1","page":"Home","title":"MapImages.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Index-1","page":"Home","title":"Index","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [MapImages]","category":"page"},{"location":"#MapImages.@latlon-Tuple{Any}","page":"Home","title":"MapImages.@latlon","text":"@allow_mapimage(funcname)\n\nmacro to make DemRsc functions equivalently defined for MapImage.\n\nThe macro just passes along the MapImage.demrsc to the function already defined \n\nUsage\n\nnearest_row(demrsc::DemRsc, lat::AbstractFloat) = Int(round(1 + (lat - d.y_first) / d.y_step))\n@allow_mapimage nearest_row\n\nThis is the same as:\n\nnearest_row(img::MapImage, lat) = nearest_row(img.demrsc, lat)\n\n\n\n\n\n","category":"macro"},{"location":"#MapImages.crop_demrsc-Tuple{Sario.DemRsc}","page":"Home","title":"MapImages.crop_demrsc","text":"crop_demrsc(demrsc::DemRsc;\n            new_rows=demrsc.rows, new_cols=demrsc.cols, \n            row_start=1, col_start=1, \n            row_step=1, col_step=1)\n\nAdjusts the old demrsc of a MapImage after slicing. Takes the 'rows' and 'cols' from demrsc and adjusts for the smaller size with a new dict\n\nThe row_start/col_start are the indices where we have sliced. This means:\n\nimg[2:3:11, :]\n\nwould have row_start = 2, row_step = 3, new_rows = 4\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.find_obj-Tuple{Base.Broadcast.Broadcasted}","page":"Home","title":"MapImages.find_obj","text":"A = find_obj(As) returns the first MapImage among the arguments.\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.get_cache_dir","page":"Home","title":"MapImages.get_cache_dir","text":"Returns the config folder for the application.  The default behavior is to return whatever is most appropriate for the operating system.\n\nThis is used to store gps timeseries data\n\nthe following folders could be returned: Mac OS X:   LibraryApplication Supportapertools Mac OS X (POSIX):   apertools Unix:   cacheapertools Unix (POSIX):   apertools\n\nforce_posix: if this is set to true then on any POSIX system the     folder will be stored in the home folder with a leading     dot instead of the XDG config home or darwin's     application support folder.\n\nBase on: https://github.com/pallets/click/blob/master/click/utils.py#L368\n\n\n\n\n\n","category":"function"},{"location":"#MapImages.grid-Tuple{Sario.DemRsc}","page":"Home","title":"MapImages.grid","text":"grid(demrsc::DemRsc; sparse=false)\n\nReturns meshgrid-like arrays X, Y\n\nExamples\n\njulia> demrsc = DemRsc(x_first=1.0, y_first=2.0, x_step=0.1, y_step=-0.2, file_length=6, width=5);\n\njulia> XX, YY = MapImages.grid(demrsc);\n\njulia> XX[1:3, 1:3]\n3×3 Array{Float64,2}:\n 1.0  1.1  1.2\n 1.0  1.1  1.2\n 1.0  1.1  1.2\n \njulia> YY[1:3, 1:3]\n3×3 Array{Float64,2}:\n 2.0  2.0  2.0\n 1.8  1.8  1.8\n 1.6  1.6  1.6\n\njulia> size(XX), size(YY)\n((6, 5), (6, 5))\n\njulia> XX, YY = MapImages.grid(demrsc, sparse=true);\n\njulia> size(XX), size(YY)\n((5,), (6,))\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.grid_extent-Tuple{Sario.DemRsc}","page":"Home","title":"MapImages.grid_extent","text":"grid_extent(demrsc::DemRsc)\n\nGet the boundaries of a grid\n\nReturns:     tuple[float]: the boundaries of the latlon grid in order:     (lonleft,lonright,latbottom,lattop)\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.intersection_corners-Tuple{Sario.DemRsc,Sario.DemRsc}","page":"Home","title":"MapImages.intersection_corners","text":"Returns:     tuple[float]: the boundaries of the intersection box of the 2 areas in order:     (lonleft,lonright,latbottom,lattop)\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.km_to_deg","page":"Home","title":"MapImages.km_to_deg","text":"Find the degrees separation for distance km assuming distance along great circle arc \n\n\n\n\n\n","category":"function"},{"location":"#MapImages.latlon_to_dist","page":"Home","title":"MapImages.latlon_to_dist","text":"Find the distance between two lat/lon points on Earth\n\nlatlon_to_dist(lat_lon_start, lat_lon_end, R=6378)\n\nUses the haversine formula: https://en.wikipedia.org/wiki/Haversine_formula so it does not account for the ellopsoidal Earth shape. Will be with about 0.5-1% of the correct value.\n\nNotes: lats and lons are in degrees, and the values used for R Earth (6373 km) are optimized for locations around 39 degrees from the equator\n\nReference: https://andrew.hedges.name/experiments/haversine/\n\n\n\n\n\n","category":"function"},{"location":"#MapImages.nsew-Tuple{Any}","page":"Home","title":"MapImages.nsew","text":"return tuple of (north, south, east, west) from demrsc data\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.read_station_llas-Tuple{}","page":"Home","title":"MapImages.read_station_llas","text":"Read in a DataFrame of gps stations with columns name, lat, lon, alt\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.rowcol_to_latlon-Tuple{Sario.DemRsc,Any,Any}","page":"Home","title":"MapImages.rowcol_to_latlon","text":"rowcol_to_latlon(demrsc::DemRsc, row, col)\n\nTakes the row, col of a pixel and finds its lat/lon within the grid\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.station_lonlat-Tuple{AbstractString}","page":"Home","title":"MapImages.station_lonlat","text":"Return the (lon, lat) of a station_name\n\n\n\n\n\n","category":"method"},{"location":"#MapImages.@allow_mapimage-Tuple{Any}","page":"Home","title":"MapImages.@allow_mapimage","text":"@allow_mapimage(funcname)\n\nmacro to make DemRsc functions equivalently defined for MapImage.\n\nThe macro just passes along the MapImage.demrsc to the function already defined \n\nUsage\n\nnearest_row(demrsc::DemRsc, lat::AbstractFloat) = Int(round(1 + (lat - d.y_first) / d.y_step))\n@allow_mapimage nearest_row\n\nThis is the same as:\n\nnearest_row(img::MapImage, lat) = nearest_row(img.demrsc, lat)\n\n\n\n\n\n","category":"macro"}]
}
