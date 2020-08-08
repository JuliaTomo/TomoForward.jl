mutable struct VolGeom{T<:AbstractFloat}
    nx::Int
    ny::Int
    nz::Int
    minX::T
    maxX::T
    minY::T
    maxY::T
    minZ::T
    maxZ::T
    spacingX::T
    spacingY::T
    spacingZ::T
end

function VolGeom(nx::Int, ny::Int, nz::Int, minX, maxX, minY, maxY, minZ, maxZ)
    spacingX = (maxX - minX) / nx
    spacingY = (maxY - minY) / ny
    spacingZ = (maxZ - minZ) / nz
    return VolGeom(nx, ny, nz, minX, maxX, minY, maxY, minZ, maxZ, spacingX, spacingY, spacingZ)
end

function VolGeom(nx::Int, ny::Int, nz::Int)
    minX = -nx / 2.0
    maxX = +nx / 2.0
    minY = -ny / 2.0
    maxY = +ny / 2.0
    minZ = -nz / 2.0
    maxZ = +nz / 2.0
    return VolGeom(nx, ny, nz, minX, maxX, minY, maxY, minZ, maxZ)
end

function VolGeom(nx::Int, ny::Int, nz::Int, minX, maxX, minY, maxY, minZ, maxZ)
    spacingX = (maxX - minX) / nx
    spacingY = (maxY - minY) / ny
    spacingZ = (maxZ - minZ) / nz
    return VolGeom(nx, ny, nz, minX, maxX, minY, maxY, minZ, maxZ, spacingX, spacingY, spacingZ)
end

"2D constructor"
function VolGeom(nx::Int, ny::Int, minX, maxX, minY, maxY)
    spacingX = (maxX - minX) / nx
    spacingY = (maxY - minY) / ny
    spacingZ = 0f0
    return VolGeom(nx, ny, 0, minX, maxX, minY, maxY, 0f0, 0f0, spacingX, spacingY, spacingZ)
end

function VolGeom(proj_geom)
    n = max(proj_geom.DetectorColCount, proj_geom.DetectorRowCount)
    half_x_width = proj_geom.DetectorColCount * proj_geom.DetectorSpacingX * 0.5
    if proj_geom.DetectorRowCount == 0
        # 2D
        vol_geom = VolGeom(n, n, -half_x_width, half_x_width, -half_x_width, half_x_width)
    else
        half_y_width = 0.5f0 * proj_geom.DetectorRowCount * proj_geom.DetectorSpacingY
        vol_geom = VolGeom(n, n, n, -half_x_width, half_x_width, -half_x_width, half_x_width, -half_y_width, half_y_width)
    end
    return vol_geom
end