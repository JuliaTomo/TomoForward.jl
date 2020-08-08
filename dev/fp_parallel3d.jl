"""
for each angle,
    compute the direction
        and choose some blocks and run CUDA code
"""

using CUDA

function fp_parallel3d_kernel()
    ny, nx, nz = size(grid)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1 # begins from 0
    jl_idx = index + 1
    hidx = Int(index % H)
    index = Int((index - hidx) / H)
    widx = Int(index % W)
    bidx = Int((index - widx) / W)
end

function fp_parallel3d(proj_geom, vol_geom)

end

nthreads = 128
totalthread = nangles * H * W
nblocks = Int(floor(totalthread / nthreads + 1)) + 1;

@cuda threads=nthreads blocks=nblocks fp_kernel!(p, surf_interp, cntcum, sz_proj, scale_bp, idx_bp)