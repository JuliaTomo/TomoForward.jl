"""
-----------------------------------------------------------------------
Copyright: 2010-2015, iMinds-Vision Lab, University of Antwerp
           2014-2015, CWI, Amsterdam
Contact: astra@uantwerpen.be
Website: http://sf.net/projects/astra-toolbox
This file is part of the ASTRA Toolbox.
The ASTRA Toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
The ASTRA Toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with the ASTRA Toolbox. If not, see <http://www.gnu.org/licenses/>.
-----------------------------------------------------------------------
"""
# https://github.com/mohamedadaly/TRex/blob/6f1b9367102c374df71d32c6ef7ede33f67c92a1/cuda/3d/par3d_fp.cu

using CUDA

const g_anglesPerBlock = 4
const g_blockSlices = 32
const g_detBlockU = 32
const g_detBlockV = 32


function fp_parallel3d_kernel!(p, vector, start_slice, start_angle, end_angle,
            nslices, ndim1, ndim2, c0, c1, c2, fx, fy, fz)
    # index = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1 # begins from 0
    # blockIdx().y

    b = c0(1,2,3)
    end_slice = nslices

    @cuprintln(b)
    # jl_idx = index + 1
    # hidx = Int(index % H)
    # index = Int((index - hidx) / H)
    # widx = Int(index % W)
    # bidx = Int((index - widx) / W)

    detX = vector[1] + vector[7] + vector[10]
    detY = vector[2] + vector[8] + vector[11]
    detZ = vector[3] + vector[9] + vector[12]

    a1 = c1(vector[4], vector[5], vector[6]) / c0(vector[4], vector[5], vector[6])
    a2 = c2(vector[4], vector[5], vector[6]) / c0(vector[4], vector[5], vector[6])
    b1 = c1(vector[4], vector[5], vector[6]) - a1 * c0(vector[4], vector[5], vector[6])
    b2 = c2(vector[4], vector[5], vector[6]) - a2 * c0(vector[4], vector[5], vector[6])

    # corr dist

    fval = 0.f0

    for detectorV=st:stend 


        # marching along dirx or diry dirz

    end

    return
end

function fp_parallel3d!(p, img, proj_geom, vol_geom)
    """
    for each angle,
        compute the direction
            and choose some blocks and run CUDA code
    """
    nangles = length(proj_geom.ProjectionAngles)
    dim_block = (g_detBlockU, g_anglesPerBlock)
    vecs = proj_geom.Vectors
    
    block_start = 1
    block_end = 1



    nx, ny, nz = vol_geom.nx, vol_geom.ny, vol_geom.nz

    f1 = (x,y,z) -> x; f2 = (x,y,z) -> y; f3 = (x,y,z) -> z

    for iangle=1:nangles+1
        theta = proj_geom.ProjectionAngles[iangle]
        dir = -1
        block_dir = 0

        # compute_raydir
        if vecs[iangle, 4] >= vecs[iangle, 5] && vecs[iangle, 4] >= vecs[iangle, 6]
            dir = 0
        elseif vecs[iangle, 5] >= vecs[iangle, 4] && vecs[iangle, 5] >= vecs[iangle, 6]
            dir = 1
        else
            dir = 2
        end

        if iangle == nangles+1 || dir != block_dir
            block_end = iangle
            if block_start != block_end
                dim_threads = (1, 1)

                if block_dir == 0
                    for i=1:g_blockSlices:nx                       
                        @cuda threads=dim_threads blocks=dim_block fp_parallel3d_kernel(p, vecs[iangle,:], i, block_start, block_end, nx, ny, nz, f1, f2, f3, f1, f2, f3)
                    end
                elseif block_dir == 1
                    c0 = (x,y,z) -> y; c1 = (x,y,z) -> x; c2 = (x,y,z) -> z
                    fx = c0; fy = c1; fz = c2;
                    for i=1:g_blockSlices:ny
                        
                        @cuda threads=dim_threads blocks=dim_block fp_parallel3d_kernel(p, vecs[iangle,:], i, block_start, block_end, ny, nx, nz, f2, f1, f3, f2, f1, f3)
                    end
                elseif block_dir == 2
                    c0 = (x,y,z) -> z; c1 = (x,y,z) -> x; c2 = (x,y,z) -> y
                    fx = (x,y,z) -> ; fy = c1; fz = c2;
                    for i=1:g_blockSlices:nz
                        @cuda threads=dim_threads blocks=dim_block fp_parallel3d_kernel(p, vecs[iangle,:], i, block_start, block_end, nz, nx, ny, f3, f1, f2, f2, f3, f1)
                    end
                end
            end            

            block_dir = dir
            block_start = iangle
        end
    end
end

nangles = 11
nthreads = 128
totalthread = nangles * 34 # * vol_geom.nx
nblocks = Int(floor(totalthread / nthreads + 1)) + 1;

fun3 = (x,y,z) -> x
# fun1(x) = x[1]

function test(A)
    return 
end

@cuda threads=(1,1) blocks=(2,2) test(proj_geom)

#------- test textures
arr = cu(zeros(3,4,5))
texture = CuTextureArray(arr)

# texarr3D = CUDA.CuTextureArray{Float32}(2,3,4)
# tex3D = CuTexture(texarr3D)
