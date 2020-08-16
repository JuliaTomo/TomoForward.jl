include("../src/proj_geom.jl")
include("../src/fp_op_common.jl")
include("../src/vol_geom.jl")
include("fp_parallel3d.jl")

using CUDA

nslice = 128
W = 128
img = zeros(Float32, W, W, nslice)
img[70:101, 40:101, 40:55] .= 1

nangles = 30
detwidth = W
spacing = 2.f0 / detwidth
proj_geom = ProjGeom(spacing, spacing, W, W, LinRange(0,pi,nangles+1)[1:nangles])

vol_geom = VolGeom(W, W, W, -1.f0, 1.f0, -1.f0, 1.f0, -1.f0, 1.f0)

p = zeros(Float32, proj_geom.DetectorRowCount, proj_geom.DetectorColCount, nangles)

img_cuda = CuArray(img)
p_cuda = CuArray(p)

@time fp_parallel3d!(p_cuda, img_cuda, proj_geom, vol_geom)

p = Array(p_cuda)
using MAT
@show maximum(p), sum(p .> 1f-8)
matwrite("../result/proj.mat", Dict("proj" => Array(p)))

