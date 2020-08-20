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
            nx, ny, nz, c0, c1, c2, detSX, detSY, detSZ, dir, scale1, scale2, out_scale)

    # project for one angle
    tid = threadIdx().x + (blockIdx.().x - 1) * blockDim().x
    I = CartesianIndices(p)
    
    
    # if tid == 1
    #     @cuprintln("$(img[0.00f0,0.0f0,0.0f0])")
    #     @cuprintln("$(img[-10.00f0,0.0f0,0.0f0])")
    #     @cuprintln("$(img[2.00f0,1.0f0,1.0f0])")
    #     @cuprintln("$(img[1.50f0,1.0f0,1.0f0])")
    # else
    #     return
    # end


    # for each detector pixel
    @inbounds if tid <= length(I)
        detV, detU = Tuple(I[tid])

        
        detX = detSX + (detU-1)*vector[7] + (detV-1)*vector[10]
        detY = detSY + (detU-1)*vector[8] + (detV-1)*vector[11]
        detZ = detSZ + (detU-1)*vector[9] + (detV-1)*vector[12]

        a1 = c1(vector[1], vector[2], vector[3]) / c0(vector[1], vector[2], vector[3])
        a2 = c2(vector[1], vector[2], vector[3]) / c0(vector[1], vector[2], vector[3])
        b1 = c1(detX, detY, detZ) - a1 * c0(detX, detY, detZ)
        b2 = c2(detX, detY, detZ) - a2 * c0(detX, detY, detZ)

        # corr dist
        corr_dist = sqrt(a1*a1*scale1 + a2*a2*scale2 + 1.f0) * out_scale

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
        # f0 -= 0.5f0; f1 -= 0.5f0; f2 -= 0.5f0; 

        # @cuprintln("$f0 $f1 $a1 $f2 $a2 $(img[f0,f1,f2])")

        value = 0.f0
        bdy = 1.f0
        for i=1:nx
            # marching along dirx or diry dirz
            if dir == 0
                if f0 >= bdy && f1 >= bdy && f0 <= nx+bdy && f1 <= ny+bdy
                    value += img[f0, f1, f2]
                end
                # value += img[f1, f0, f2]
            elseif dir == 1
                if f0 >= bdy && f1 >= bdy && f1 <= nx+bdy && f0 <= ny+bdy
                    value += img[f1, f0, f2]
                end
                # value += img[f0, f1, f2]
            else
                # in parallel beam, typically dir \neq 2
                # TODO
                # value += img[f1, f2, f0]
            end
            f0 += 1.0f0;
			f1 += a1;
            f2 += a2;
            
            # debugging
            if img[f0,f1,f2] > 0.f0
                # @cuprintln(img[127.f0, 128.f0, 127.f0], " ", img[127.9f0, 128.f0, 127.f0], " ", img[0.4f0, 0.5f0, 0.5f0], " ", img[-1.5f0, -1.5f0, -1.5f0])
                # @cuprintln("$f0, $f1, $f2, $(img[f0,f1,f2])")
            end
        end

        value *= corr_dist
        p[detV, detU] = value
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
    elseif size(img, 1) != vol_geom.ny || size(img, 2) != vol_geom.nx || size(img, 3) != vol_geom.nz
        # error("mismatch between img size and volume geometry")
    end

    nthreads = 512
    nblocks = cld(size(p,1)*size(p,2), nthreads)

    nangles = length(proj_geom.ProjectionAngles)
    vecs = proj_geom.Vectors
    
    nx, ny, nz = vol_geom.nx, vol_geom.ny, vol_geom.nz
    # out_scale = 1.f0

    texturearray = CuTextureArray(img)
    # Create a texture object and bind it to the texture memory created above
    img_texture = CuTexture(texturearray; filter_mode=CUDA.CU_TR_FILTER_MODE_LINEAR)
    nrow = proj_geom.DetectorRowCount
    ncol = proj_geom.DetectorColCount

    for iangle=1:nangles
        out_scale = vol_geom.spacingX

        vector = vecs[iangle,:]
        # translation https://github.com/astra-toolbox/astra-toolbox/blob/10d87f45bc9311c0408e4cacdec587eff8bc37f8/include/astra/GeometryUtil3D.h
        
        # scaling
        vector[[1,7,10]] ./= vol_geom.spacingX
        vector[[2,8,11]] ./= vol_geom.spacingY
        vector[[3,9,12]] ./= vol_geom.spacingZ
        
        detSX = vector[4] -0.5f0*vector[7]*ncol -0.5f0*vector[10]*nrow
        detSY = vector[5] -0.5f0*vector[8]*ncol -0.5f0*vector[11]*nrow
        detSZ = vector[6] -0.5f0*vector[9]*ncol -0.5f0*vector[12]*nrow

        if false
            # translation for the non-center volume
            dx = -(vol_geom.minX + vol_geom.maxX) / 2.f0;
            dy = -(vol_geom.minX + vol_geom.maxX) / 2.f0;
            dz = -(vol_geom.minX + vol_geom.maxX) / 2.f0;

            detSX += dx; detSY += dy; detSZ += dz;
        end

        # detSX /= vol_geom.spacingX
        # detSY /= vol_geom.spacingY
        # detSZ /= vol_geom.spacingZ

        detSX += 0.5f0*vector[7] + 0.5f0*vector[10]
        detSY += 0.5f0*vector[8] + 0.5f0*vector[11]
        detSZ += 0.5f0*vector[9] + 0.5f0*vector[12]   
    
        p_ = view(p, :, :, iangle)
        # compute_raydir
        if abs(vector[1]) >= abs(vector[2]) && abs(vector[1]) >= abs(vector[3])
            dir = 0
            scale1 = (vol_geom.spacingY / vol_geom.spacingX)^2
            scale2 = (vol_geom.spacingZ / vol_geom.spacingX)^2
        elseif abs(vector[2]) >= abs(vector[1]) && abs(vector[2]) >= abs(vector[3])
            dir = 1
            scale1 = (vol_geom.spacingX / vol_geom.spacingY)^2
            scale2 = (vol_geom.spacingZ / vol_geom.spacingY)^2
            out_scale = vol_geom.spacingY
        else
            # dir = 2
            # scale1 = (vol_geom.spacingX / vol_geom.spacingZ)^2
            # scale2 = (vol_geom.spacingY / vol_geom.spacingZ)^2
        end

        vector = cu(vector)

        if dir == 0
            c0 = (x,y,z) -> x; c1 = (x,y,z) -> y; c2 = (x,y,z) -> z
            @cuda threads=nthreads blocks=nblocks fp_parallel3d_kernel!(p_, img_texture, vector, iangle, nx, ny, nz, c0, c1, c2, detSX, detSY, detSZ, dir, scale1, scale2, out_scale)
        elseif dir == 1
            c0 = (x,y,z) -> y; c1 = (x,y,z) -> x; c2 = (x,y,z) -> z
            @cuda threads=nthreads blocks=nblocks fp_parallel3d_kernel!(p_, img_texture, vector, iangle, ny, nx, nz, c0, c1, c2, detSX, detSY, detSZ, dir, scale1, scale2, out_scale)
        elseif dir == 2
            # c0 = (x,y,z) -> z; c1 = (x,y,z) -> x; c2 = (x,y,z) -> y
            # @cuda threads=nthreads blocks=nblocks fp_parallel3d_kernel!(p_, img_texture, vector, iangle, nz, nx, ny, c0, c1, c2, detSX, detSY, detSZ, dir, scale1, scale2, out_scale)
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