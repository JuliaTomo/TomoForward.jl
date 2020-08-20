ENV["PYTHON"] = raw"/nobackup/jakoo/miniconda3/envs/ctdr/bin/python"
using PyCall
np = pyimport("numpy")
astra = pyimport("astra")
@show astra.__version__
nslice = 128
W = 100

nangles = 30
angles = Float32.(LinRange(0,pi,nangles+1)[1:nangles])
detwidth = W
# spacing = 1.f0
spacing = 2.f0 / detwidth
pg = astra.create_proj_geom("parallel3d", spacing, spacing, W, W, angles)

if spacing < 1.f0
    vg = astra.create_vol_geom(W, W, nslice, -1.f0, 1.f0, -1, 1, -1, 1)
else
    vg = astra.create_vol_geom(W, W, nslice)
end
    # Coordinate order: slice, row, column (z, y, x)
img_np = np.zeros((nslice, W, W), dtype=np.float32)
img_np[40:W, 40:W-10, 70:W-10] .= 1
img_np[20, 2:11, 13:40] .= 3
# img_np[20:40, 2:11, 3:4] .= 2
# img_np[20:40, :, :] .= 2
# img_np .= 1
# img_np[64,64,64] = 1


@time pid, proj_ = astra.create_sino3d_gpu(img_np, pg, vg)

proj = permutedims(proj_, (1,3,2))
@show maximum(proj), sum(proj .> 1f-8)
# proj = mapslices(x -> reshape(vec(x),W,W), proj, dims=[1,2])

if true
    # s0 = astra.data3d.create("-proj3d", pg)
    # Initialized to a matrix:
    # Coordinate order: row (v), angle, column (u)
    # A = np.ones((W, nangles, W))
    A = proj_
    pid = astra.data3d.create("-proj3d", pg, A)
end

@time iid, img_bp_astra = astra.create_backprojection3d_gpu(pid, pg, vg)
@show maximum(img_bp_astra)
img_bp_astra_ = permutedims(img_bp_astra, (3,2,1))
using MAT
matwrite("../result/proj_astra.mat", Dict("proj_astra" => proj, "img_bp_astra" => img_bp_astra_))

