# module proj_geom

include("util_astra.jl")
export ProjGeom

"""Projection geometry structure

    Type::String
    DetectorSpacingX::AbstractFloat
    DetectorSpacingY::AbstractFloat
    DetectorColCount::Int
    DetectorRowCount::Int
    ProjectionAngles::Array{AbstractFloat,1}
    Vectors::Array{AbstractFloat,2}

    # parallel 2d
    ProjGeom(spacing, detcount, angles) = new("parallel2d", spacing, 0.0, detcount, 0, angles, geom_2vec_parallel2d(b, c))
    
    # parallel 3d
    ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles) = new("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, geom_2vec_parallel3d(detspacingx, detspacingy, angles))

    # cone
    ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles, source_origin, origin_det) = new("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, geom_2vec_cone(detspacingx, detspacingy, angles, source_origin, origin_det))
"""
struct ProjGeom
    Type::String
    DetectorSpacingX::AbstractFloat # (@don't need to save) already inherited in Vectors
    DetectorSpacingY::AbstractFloat
    DetectorColCount::Int
    DetectorRowCount::Int
    ProjectionAngles::Array{AbstractFloat,1}
    Vectors::Array{AbstractFloat,2}
        
    # parallel 2d
    ProjGeom(spacing, detcount, angles) = new("parallel2d", spacing, 0.0, detcount, 0, angles, geom_2vec_parallel2d(spacing, angles))
    
    # parallel 3d
    ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles) = new("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, angles, geom_2vec_parallel3d(detspacingx, detspacingy, angles))
    
    # cone
    ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles, source_origin, origin_det) = new("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, angles, geom_2vec_cone(detspacingx, detspacingy, angles, source_origin, origin_det))
    
end

# end