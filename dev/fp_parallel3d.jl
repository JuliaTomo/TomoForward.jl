"""

"""

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
#https://github.com/astra-toolbox/astra-toolbox/blob/06322245a638e435d29bc8f027aed0bd976139d0/cuda/3d/par3d_fp.cu

using CUDA

function fp_parallel3d_kernel!(p, img, vector, iangle,
            nx, ny, nz, c0, c1, c2, detSX, detSY, detSZ, dir, scale)

    # project for one angle
    tid = threadIdx().x + (blockIdx.().x - 1) * blockDim().x
    I = CartesianIndices(p)
    
    # for each detector pixel
    @inbounds if tid <= length(I)
        detV, detU = Tuple(I[tid])
        
        detX = detSX + detU*vector[7] + detV*vector[10]
        detY = detSY + detU*vector[8] + detV*vector[11]
        detZ = detSZ + detU*vector[9] + detV*vector[12]

        a1 = c1(vector[1], vector[2], vector[3]) / c0(vector[1], vector[2], vector[3])
        a2 = c2(vector[1], vector[2], vector[3]) / c0(vector[1], vector[2], vector[3])
        b1 = c1(detX, detY, detZ) - a1 * c0(detX, detY, detZ)
        b2 = c2(detX, detY, detZ) - a2 * c0(detX, detY, detZ)

        # corr dist
        corr_dist = sqrt(a1*a1 + a2*a2 + 1.f0) * scale
        value = 0.f0

        """
        x =  1 x +  0
        y = a1 x + b1
        z = a2 x + b2
        """
        # start_slice = 0

        f0 = + 0.5f0
        f1 = a1 * ( - 0.5f0*nx + 0.5f0) + b1 # y in world coordinate
        f1 += 0.5f0*ny - 0.5f0 + 0.5f0 # y in raster coordinate
        f2 = a2 * ( - 0.5f0*nx + 0.5f0) + b2
        f2 += 0.5f0*nz - 0.5f0 + 0.5f0

        f0 += 1.f0; f1 += 1.f0; f2 += 1.f0; # Julia indexing

        # @cuprintln("$f0 $f1 $a1 $f2 $a2 $(img[f0,f1,f2])")

        for i=1:nx
            # marching along dirx or diry dirz
            if dir == 0
                value += img[f0, f1, f2]
            elseif dir == 1
                value += img[f1, f0, f2]
            else
                value += img[f1, f2, f0]
            end
            f0 += 1.0f0;
			f1 += a1;
			f2 += a2;
        end

        value *= corr_dist
        p[tid] = value
    end
    return
end

"""
    fp_parallel3d!(p, img, proj_geom, vol_geom)

forward projection of 3D volume.
"""
function fp_parallel3d!(p, img, proj_geom, vol_geom)
    """
    for each angle,
        compute the direction
            and choose some blocks and run CUDA code
    """

    if abs(vol_geom.minX + vol_geom.maxX) > 1e-8
        error("not implemented yet for non-cube volume")
    end

    nthreads = 128
    nblocks = cld(size(p,1)*size(p,2), nthreads)

    nangles = length(proj_geom.ProjectionAngles)
    vecs = cu(proj_geom.Vectors)
    
    nx, ny, nz = vol_geom.nx, vol_geom.ny, vol_geom.nz

    texturearray = CuTextureArray(img)
    # Create a texture object and bind it to the texture memory created above
    img_texture = CuTexture(texturearray)
    nrow = proj_geom.DetectorRowCount
    ncol = proj_geom.DetectorColCount

    for iangle=1:nangles
        vector = vecs[iangle,:]
        vector[[1,7,10]] ./= vol_geom.spacingX
        vector[[2,8,11]] ./= vol_geom.spacingY
        vector[[3,9,12]] ./= vol_geom.spacingZ
        scale = vol_geom.spacingX

        detSX = vector[4] - 0.5f0*vector[7]*ncol - 0.5f0*vector[10]*nrow
        detSY = vector[5] -0.5f0*vector[8]*ncol -0.5f0*vector[11]*nrow
        detSZ = vector[6] -0.5f0*vector[9]*ncol -0.5f0*nrow*vector[12]

        detSX += 0.5f0*vector[7] + 0.5f0*vecs[10]
        detSY += 0.5f0*vecs[8] + 0.5f0*vecs[11]
        detSZ += 0.5f0*vecs[9] + 0.5f0*vecs[12]   
    
        p_ = view(p, :, :, iangle)
        # compute_raydir
        if abs(vector[1]) >= abs(vector[2]) && abs(vector[1]) >= abs(vector[3])
            dir = 0
        elseif abs(vector[2]) >= abs(vector[1]) && abs(vector[2]) >= abs(vector[3])
            dir = 1
        else
            dir = 2
        end

        if dir == 0
            c0 = (x,y,z) -> x; c1 = (x,y,z) -> y; c2 = (x,y,z) -> z
            # fx = c0; fy = c1; fz = c2;
            @cuda threads=nthreads blocks=nblocks fp_parallel3d_kernel!(p_, img_texture, vector, iangle, nx, ny, nz, c0, c1, c2, detSX, detSY, detSZ, dir, scale)
        elseif dir == 1
            c0 = (x,y,z) -> y; c1 = (x,y,z) -> x; c2 = (x,y,z) -> z
            @cuda threads=nthreads blocks=nblocks fp_parallel3d_kernel!(p_, img_texture, vector, iangle, ny, nx, nz, c0, c1, c2, detSX, detSY, detSZ, dir, scale)
        elseif dir == 2
            c0 = (x,y,z) -> z; c1 = (x,y,z) -> x; c2 = (x,y,z) -> y
            @cuda threads=nthreads blocks=nblocks fp_parallel3d_kernel!(p_, img_texture, vector, iangle, nz, nx, ny, c0, c1, c2, detSX, detSY, detSZ, dir, scale)
        end
    end
end

# if iangle == nangles+1 || dir != block_dir
        #     block_end = iangle
        #     if block_start != block_end
        #         dim_threads = (1, 1)

        #         if block_dir == 0
        #             for i=1:g_blockSlices:nx                       
        #                 @cuda threads=dim_threads blocks=dim_block fp_parallel3d_kernel!(p, img, vecs[iangle,:], i, block_start, block_end, nx, ny, nz, f1, f2, f3, f1, f2, f3)
        #             end
        #         elseif block_dir == 1
        #             c0 = (x,y,z) -> y; c1 = (x,y,z) -> x; c2 = (x,y,z) -> z
        #             fx = c0; fy = c1; fz = c2;
        #             for i=1:g_blockSlices:ny
                        
        #                 @cuda threads=dim_threads blocks=dim_block fp_parallel3d_kernel!(p, img, vecs[iangle,:], i, block_start, block_end, ny, nx, nz, f2, f1, f3, f2, f1, f3)
        #             end
        #         elseif block_dir == 2
        #             c0 = (x,y,z) -> z; c1 = (x,y,z) -> x; c2 = (x,y,z) -> y
        #             fx = (x,y,z) -> ; fy = c1; fz = c2;
        #             for i=1:g_blockSlices:nz
        #                 @cuda threads=dim_threads blocks=dim_block fp_parallel3d_kernel!(p, img, vecs[iangle,:], i, block_start, block_end, nz, nx, ny, f3, f1, f2, f2, f3, f1)
        #             end
        #         end
        #     end            

        #     block_dir = dir
        #     block_start = iangle