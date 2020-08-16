using PyCall
np = pyimport("numpy")
astra = pyimport("astra")
nslice = 128
W = 128

nangles = 30
angles = Float32.(LinRange(0,pi,nangles+1)[1:nangles])
detwidth = W
spacing = 2.f0 / detwidth
pg = astra.create_proj_geom("parallel3d", spacing, spacing, W, W, angles)
vg = astra.create_vol_geom(W, W, W, -1.f0, 1.f0, -1, 1, -1, 1)

# Coordinate order: slice, row, column (z, y, x)
img_np = np.zeros((W, W, nslice), dtype=np.float32)
img_np[40:55, 40:101, 70:101] .= 1


@time pid, proj_ = astra.create_sino3d_gpu(img_np, pg, vg)

proj = permutedims(proj_, (1,3,2))
@show maximum(proj), sum(proj .> 1f-8)
# proj = mapslices(x -> reshape(vec(x),W,W), proj, dims=[1,2])

using MAT
matwrite("../result/proj_astra.mat", Dict("proj_astra" => proj))
