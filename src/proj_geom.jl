include("util_astra.jl")

"""Projection geometry structure

    struct ProjGeom {T<:AbstractFloat}
    Type::String
    DetectorSpacingX::T 
    DetectorSpacingY::T
    DetectorColCount::Int
    DetectorRowCount::Int
    ProjectionAngles::Array{T,1}
    Vectors::Array{T,2}

    # parallel 2d
    ProjGeom(spacing, detcount, angles) = new("parallel2d", spacing, 0.0, detcount, 0, angles, geom_2vec_parallel2d(b, c))
    
    # parallel 3d
    ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles) = new("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, geom_2vec_parallel3d(detspacingx, detspacingy, angles))

    # cone
    ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles, source_origin, origin_det) = new("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, geom_2vec_cone(detspacingx, detspacingy, angles, source_origin, origin_det))
"""
struct ProjGeom{T<:AbstractFloat}
    Type::String
    DetectorSpacingX::T 
    DetectorSpacingY::T
    DetectorColCount::Int
    DetectorRowCount::Int
    ProjectionAngles::Array{T,1}
    Vectors::Array{T,2}
end

# parallel 2d
ProjGeom(spacing, detcount, angles) = new("parallel2d", spacing, 0.0, detcount, 0, angles, geom_2vec_parallel2d(spacing, angles))
    
# parallel 3d
ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles) = new("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, angles, geom_2vec_parallel3d(detspacingx, detspacingy, angles))

# cone
ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles, source_origin, origin_det) = new("cone", detspacingx, detspacingy, detcolcount, detrowcount, angles, geom_2vec_cone(detspacingx, detspacingy, angles, source_origin, origin_det))
