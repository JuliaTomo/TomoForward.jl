using TomoForward

N = 128
img = zeros(N, N)
img[70:101, 40:101] .= 1

nangles = 200
detcount = N
# detcount = Int(floor(size(img,1)*1.4))
deetcount = 512
src_origin = 1500.0
det_origin = 500.0
detspacing = (src_origin+det_origin) / src_origin
angles = LinRange(0,pi,nangles+1)[1:nangles]
proj_geom = ProjGeomFan(detspacing, detcount, angles, src_origin, det_origin)

# test line projection model
A = fp_op_fan_line(proj_geom, size(img, 1), size(img, 2))
p = A * vec(img);
bp = A' * vec(p)
bp_img = reshape(bp, size(img, 1), size(img, 2))

using PyPlot
title("back projection with line model")
imshow(bp_img)