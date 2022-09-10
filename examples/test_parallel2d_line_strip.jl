using TomoForward

img = zeros(200, 200)
img[70:101, 40:101] .= 1

nangles = 200
detcount = Int(floor(size(img,1)*1.4))
proj_geom = ProjGeom(1.0, detcount, LinRange(0,pi,nangles+1)[1:nangles])

# test line projection model
A = fp_op_parallel2d_line(proj_geom, size(img, 1), size(img, 2))
p = A * vec(img);
bp = A' * vec(p)
bp_img = reshape(bp, size(img, 1), size(img, 2))

using PyPlot
title("back projection with line model")
imshow(bp_img)

# test strip projection model
A_strip = fp_op_parallel2d_strip(proj_geom, size(img, 1), size(img, 2))
ps = A_strip * vec(img);
bps = A_strip' * vec(p)
bps_img = reshape(bps, size(img, 1), size(img, 2))

title("back projection with strip model")
imshow(bps_img)