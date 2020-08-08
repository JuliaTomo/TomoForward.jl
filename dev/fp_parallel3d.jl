"""
for each angle,
    compute the direction
        and choose some blocks and run CUDA code
"""

using CUDA

function fp_parallel3d_kernel()
end

function fp_parallel3d(inters_pos, inters_idx, inters_dir, inters_t, cnt_collide, cntcum, grid, orig, dir,is_first_pass::Int, H::Int, W::Int, nangles::Int,
    minX::Float32, maxX::Float32, minY::Float32, maxY::Float32, minZ::Float32, maxZ::Float32,
    spacingX::Float32, spacingY::Float32, spacingZ::Float32)

    ny, nx, nz = size(grid)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1 # begins from 0
    jl_idx = index + 1
    hidx = Int(index % H)
    index = Int((index - hidx) / H)
    widx = Int(index % W)
    bidx = Int((index - widx) / W)

end

nthreads = 128
totalthread = nangles * H * W
nblocks = Int(floor(totalthread / nthreads + 1)) + 1;

@cuda threads=nthreads blocks=nblocks fp_kernel!(p, surf_interp, cntcum, sz_proj, scale_bp, idx_bp)