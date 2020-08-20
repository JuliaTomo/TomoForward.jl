include("../src/vol_geom.jl")
include("../src/proj_geom.jl")
include("../src/fp_op_common.jl")
include("fp_parallel3d.jl")
include("bp_parallel3d.jl")


using CUDA

nslice = 128
W = 100
img = zeros(Float32, W, W, nslice)
img[70:W-10, 40:W-10, 40:W] .= 1
img[13:40, 2:11, 20] .= 3
# img = ones(Float32, W, W, nslice)
# img[1,1,1] = 10.f0
# img[2,1,1] = 20
# img[64,64,64] = 1

nangles = 30
detwidth = W
# spacing = 1.f0
spacing = 2.f0 / detwidth
# 
proj_geom = ProjGeom(spacing, spacing, W, W, LinRange(0,pi,nangles+1)[1:nangles])

if spacing < 1.f0
    vol_geom = VolGeom(W, W, nslice, -1.f0, 1.f0, -1.f0, 1.f0, -1.f0, 1.f0)
else
    @show spacing
    vol_geom = VolGeom(W, W, nslice)
end
p = zeros(Float32, proj_geom.DetectorRowCount, proj_geom.DetectorColCount, nangles)

img_cuda = CuArray(img)
p_cuda = CuArray(p)

@time fp_parallel3d!(p_cuda, img_cuda, proj_geom, vol_geom)
p = Array(p_cuda)
@show maximum(p), sum(p .> 1f-8)

if true
    # p_cuda .= 1.0f0
    # p_cuda[10:50,20:30,10:11] .= 1.f0
    # p_cuda[10:30,:,:] = 1.0
end

img_bp = cu(zeros(size(img_cuda)))

@time bp_parallel3d!(img_bp, p_cuda, proj_geom, vol_geom)

@show maximum(img_bp)

using MAT
matwrite("../result/proj.mat", Dict("proj" => p, "img_bp" => Array(img_bp)))


@show sum(abs2.(p .- proj)) / length(p)
@show sum(abs2.(Array(img_bp) .- img_bp_astra_)) / length(img_bp_astra_)

