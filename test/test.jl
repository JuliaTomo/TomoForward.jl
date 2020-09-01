function test_astra()
    astra = pyimport("astra")
    vg = astra.create_vol_geom(vol_geom.ny, vol_geom.nx, -1, 1, -1, 1)
    pg = astra.create_proj_geom("parallel", proj_geom.DetectorSpacingX, proj_geom.DetectorColCount, proj_geom.ProjectionAngles)
    proj_id = astra.create_projector("line", pg, vg)
    sinogram_id, sinogram = astra.create_sino(img_gt, proj_id)
    # sinogram = reshape(vec(sinogram), size(sinogram, 2), :)
    # sinogram = permutedims(sinogram, (2, 1))
    p_astra = vec(sinogram)

end
