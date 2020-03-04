using Revise
import TomoForward
# include("../src/simulation/build_octtree.jl")

using FileIO
using MeshIO

mesh = load("../test/test_data/bunny.obj")

vs_ = mesh.vertices
min_ = minimum(vs_)
max_ = maximum(vs_)
vs = map(x -> (x - min_) ./ (max_ - min_) .* 2 .- 1 , vs_)

H, W = 100, 100
# H, W = 50, 50 # detector size
nangles = 1
angles = LinRange(0, pi, nangles+1)[1:end-1]
proj_geom = TomoForward.ProjGeom(1.9/W, 1.9/H, H, W, angles)

sinogram, tt = TomoForward.fp_mesh(proj_geom, vs, mesh.faces)

# include("../src/simulation/util_octtree.jl")
# include("../src/simulation/raytracing.jl")
# include("../src/simulation/fp_mesh.jl")
# src_pos = [0,0.,0]
# ray_dir=[0,1.,0]
# orig = SVector(1.0, 1.0, 1.0) * -1.0
# root = Cell(orig, SVector(2.0, 2, 2), [])
# build_octtree(root, vs, mesh.faces)

# faces_candidate = get_face_candidate(src_pos, ray_dir, root)
using PyPlot
idx=1
imshow(sinogram[idx,:,:])
show()