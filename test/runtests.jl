using MapImages
import Sario: DemRsc
using Test
using Documenter

@testset "MapImages.jl" begin
    doctest(MapImages)


    test_data = reshape(collect(1:30), (6, 5))
    demrsc = DemRsc(x_first=1.0, y_first=2.0, x_step=0.1, y_step=-0.2, file_length=6, width=5)
    mapim = MapImage(test_data, demrsc)
    out = MapImages.crop_demrsc(demrsc, new_rows=2, new_cols=2)

    @test out.width == 2
    @test out.file_length == 2


    out = MapImages.crop_demrsc(demrsc, row_start=2, col_start=2, new_rows=2, new_cols=2)
    @test out.x_first == 1.1
    @test out.y_first == 1.8

    out = MapImages.crop_demrsc(demrsc, new_rows=2, new_cols=2, row_step=2, col_step=2)
    @test out.x_step == 0.2
    @test out.y_step == -0.4


end
