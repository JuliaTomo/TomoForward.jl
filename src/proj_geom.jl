include("util_astra.jl")

"""Projection geometry structure

    struct ProjGeom{T<:AbstractFloat}
        Type::String
        DetectorSpacingX::T 
        DetectorSpacingY::T
        DetectorColCount::Int
        DetectorRowCount::Int
        ProjectionAngles::Array{T,1}
        Vectors::Array{T,2}
        DistanceOriginSource
        DistanceOriginDetector
    
    "parallel 2d"
    ProjGeom(spacing, detcount, angles) = ProjGeom("parallel2d", spacing, 0.0, detcount, 0, Array{Float64}(angles), geom_2vec_parallel2d(spacing, angles), nothing, nothing)

    "fan"
    ProjGeom(spacing, detcount, angles, src_origin, det_origin) = ProjGeom("fan", spacing, 0.0, detcount, 0, Array{Float64}(angles), geom_2vec_fan(spacing, angles, src_origin, det_origin), src_origin, det_origin)

    "parallel 3d"
    ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles) = ProjGeom("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, Array{Float64}(angles), geom_2vec_parallel3d(detspacingx, detspacingy, angles), nothing, nothing)
    
    "cone"
    ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles, source_origin, origin_det) = ProjGeom("cone", detspacingx, detspacingy, detcolcount, detrowcount, Array{Float64}(angles), geom_2vec_cone(detspacingx, detspacingy, angles, source_origin, origin_det), source_origin, origin_det)
"""
mutable struct ProjGeom{T<:AbstractFloat}
    Type::String
    DetectorSpacingX::T 
    DetectorSpacingY::T
    DetectorColCount::Int
    DetectorRowCount::Int
    ProjectionAngles::Array{T,1}
    Vectors::Array{T,2}
end

"parallel 2d"
ProjGeom(spacing, detcount, angles) = ProjGeom("parallel2d", spacing, 0.0, detcount, 0, Array(angles), geom_2vec_parallel2d(spacing, angles))

"parallel 3d"
ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles) = ProjGeom("parallel3d", detspacingx, detspacingy, detcolcount, detrowcount, Array(angles), geom_2vec_parallel3d(detspacingx, detspacingy, angles))

"fan"
ProjGeomFan(spacing, detcount, angles, src_origin, det_origin) = ProjGeom("fan", spacing, 0.0, detcount, 0, Array(angles), geom_2vec_fan(spacing, angles, src_origin, det_origin))

"cone"
ProjGeom(detspacingx, detspacingy, detrowcount, detcolcount, angles, source_origin, origin_det) = ProjGeom("cone", detspacingx, detspacingy, detcolcount, detrowcount, Array(angles), geom_2vec_cone(detspacingx, detspacingy, angles, source_origin, origin_det))