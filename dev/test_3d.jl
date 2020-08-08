include("../src/proj_geom.jl")
include("../src/fp_op_common.jl")

using CUDA

nslice = 128
img = ones(Float32, 128, 128, nslice)
#img[70:101, 40:101, 40:55] .= 1

nangles = 10
proj_geom_ = ProjGeom(1.0, 1.0, nslice, 128, LinRange(0,pi,nangles+1)[1:nangles])

p = zeros(Float32, proj_geom_.DetectorRowCount, proj_geom_.DetectorColCount, nangles)

img_cuda = CuArray(img)
p_cuda = CuArray(p)

nthreads = 128 # 1024
nblocks = div(nangles*proj_geom_.DetectorRowCount, nthreads) + 1

println("nthreads $nthreads, nblocks $nblocks")

Vectors = CuArray(Array{Float32,2}(proj_geom_.Vectors))

H, W, Z = size(img)

maxX=W
minX=0
minY=H
maxY=0

pixelspacingX = (maxX-minX) / W
pixelspacingY = (maxY-minY) / H

invpixelspacingX = 1. / pixelspacingX
invpixelspacingY = 1. / pixelspacingY

Ex = minX + pixelspacingX*0.5 # min x
Ey = maxY - pixelspacingY*0.5 # max y

pixelspacingX = Float32(pixelspacingX)
pixelspacingY = Float32(pixelspacingY)

Ex = Float32(Ex)
Ey = Float32(Ey)

@cuda blocks=nblocks threads=nthreads fp_parallel3d_strip_kernel(p_cuda, img_cuda, Vectors, nangles, proj_geom_.DetectorRowCount, proj_geom_.DetectorColCount, Ex, Ey, pixelspacingX, pixelspacingY)

