module TomoForward

# common
include("proj_geom.jl")
export ProjGeom, ProjGeomFan

include("vol_geom.jl")
export VolGeom

# FP
include("fp_op_common.jl")
include("fp_op_parallel2d_line.jl")
include("fp_op_parallel2d_strip.jl")
include("fp_op_fan_line.jl")
export fp_op_parallel2d_line, fp_op_parallel2d_strip, fp_op_fan_line


include("fp_op_utils.jl")
export fp_2d_slices!

# util_astra
export geom_2vec_parallel2d, geom_2vec_parallel3d, geom_2vec_fan, geom_2vec_cone

end
