using TomoForward

img = zeros(200, 200)
img[70:101, 40:101] .= 1

nangles = 200
detcount = Int(floor(size(img,1)*1.4))
src_origin = 800.0
det_origin = 500.0
angles = LinRange(0,pi,nangles+1)[1:nangles]
proj_geom = ProjGeom(1.0, detcount, Array(angles), src_origin, det_origin)

# test line projection model
A = fp_op_fan_line(proj_geom, size(img, 1), size(img, 2))
p = A * vec(img);
bp = A' * vec(p)
bp_img = reshape(bp, size(img, 1), size(img, 2))

using PyPlot
title("back projection with line model")
imshow(bp_img)