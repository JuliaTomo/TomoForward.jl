using LinearAlgebra

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
# https://github.com/astra-toolbox/astra-toolbox/blob/master/cuda/3d/par3d_bp.cu

using CUDA

"""
kernel function 
thread : x and y index for the volume

# Args
- cons : constant values
"""
function bp_parallel3d_kernel!(img, p, H, W, nangles, nx, ny, nz, cc, scale)

    # project for one angle
    tid = threadIdx().x + (blockIdx.().x - 1) * blockDim().x
    I = CartesianIndices((1:ny, 1:nz))
    
    # for each detector pixel
    @inbounds if tid <= length(I)
        bdy = 1f0
        iy, iz = Tuple(I[tid])
        ix = 1

        # xyz in world coordinate
        fx = ix-1 - 0.5f0*nx + 0.5f0
        fy = iy-1 - 0.5f0*ny + 0.5f0
        fz = iz-1 - 0.5f0*nz + 0.5f0

        for iangle = 1:nangles
            fu = cc[4,iangle] + fx*cc[1,iangle] + fy*cc[2,iangle] + fz*cc[3,iangle]
            fv = cc[8,iangle] + fx*cc[5,iangle] + fy*cc[6,iangle] + fz*cc[7,iangle]
            fangle = iangle + 0.5f0

            # fv -= 0.5f0; fu -= 0.5f0;
            fv += 1.f0; fu += 1.f0; # julia indexing
            
            for ix=1:nx
                if fu >= bdy && fv >= bdy && fu <= W+bdy && fv <= H+bdy
                    img[ix, iy, iz] += p[fv, fu, fangle] * cc[9,iangle] * scale
                end

                fu += cc[1,iangle]
                fv += cc[5,iangle]
                if iangle == 1
                    # @cuprintln(fu, " ", fv)
                end
            end
            
        end
    end
    return
end

det3x(b,c) =  b[2]*c[3] - b[3]*c[2]
det3y(b,c) = -(b[1]*c[3] - b[3]*c[1])
det3z(b,c) =  b[1]*c[2] - b[2]*c[1]
det3(a,b,c) = a[1]*det3x(b,c)+a[2]*det3y(b,c)+a[3]*det3z(b,c)
cross3(a,b) = Array{Float32}(det3x(a,b), det3y(a,b), det3z(a,b))

function scaled_cross3(a, b, sc)
    ret = cross(a,b)
    ret[1] *= sc[2]*sc[3]
    ret[2] *= sc[1]*sc[3]
    ret[3] *= sc[1]*sc[2]
    return ret
end

function compute_constant(proj_geom, vol_geom)
    
    nangles = length(proj_geom.ProjectionAngles)
    vecs = proj_geom.Vectors

    ncol = proj_geom.DetectorColCount
    nrow = proj_geom.DetectorRowCount
    
    out = zeros(Float32, 9, nangles)
    for i=1:nangles
        # https://github.com/astra-toolbox/astra-toolbox/blob/10d87f45bc9311c0408e4cacdec587eff8bc37f8/include/astra/GeometryUtil2D.h
        vector = proj_geom.Vectors[i,:]

        detSX = vector[4] -0.5f0*vector[7]*ncol -0.5f0*vector[10]*nrow
        detSY = vector[5] -0.5f0*vector[8]*ncol -0.5f0*vector[11]*nrow
        detSZ = vector[6] -0.5f0*vector[9]*ncol -0.5f0*vector[12]*nrow

        vector[[1,7,10]] ./= vol_geom.spacingX
        vector[[2,8,11]] ./= vol_geom.spacingY
        vector[[3,9,12]] ./= vol_geom.spacingZ
        
        detSX /= vol_geom.spacingX
        detSY /= vol_geom.spacingY
        detSZ /= vol_geom.spacingZ

        u = vector[7:9]
        v = vector[10:12]
        r = vector[1:3]
        d = [detSX, detSY, detSZ]

        den = det3(r,u,v)
        out[1,i] = -det3x(r,v) / den
        out[2,i] = -det3y(r,v) / den
        out[3,i] = -det3z(r,v) / den
        out[4,i] = -det3(r,d,v) / den
        out[5,i] = det3x(r,u) / den
        out[6,i] = det3y(r,u) / den
        out[7,i] = det3z(r,u) / den
        out[8,i] = det3(r,d,u) / den
        
        # s
        out[9,i] = 1.f0 / norm(scaled_cross3(u,v,[vol_geom.spacingX, vol_geom.spacingY, vol_geom.spacingZ]))
    end
    return out
end

"""
    fp_parallel3d!(p, img, proj_geom, vol_geom)

forward projection of 3D volume.
"""
function bp_parallel3d!(img, p, proj_geom, vol_geom)
    """
    for each angle,
        compute the direction
            and choose some blocks and run CUDA code
    """
    global cons
    # if size(img, 1) != vol_geom.ny || size(img, 2) != vol_geom.nx || size(img, 3) != vol_geom.nz
        # error("mismatch between img size and volume geometry")
    # end
    scale = vol_geom.spacingX*vol_geom.spacingY*vol_geom.spacingZ

    nthreads = 1024
    nblocks = cld(size(img,2)*size(img,3), nthreads)

    cons = compute_constant(proj_geom, vol_geom)
    cons = cu(cons)

    nangles = length(proj_geom.ProjectionAngles)
    nx, ny, nz = vol_geom.nx, vol_geom.ny, vol_geom.nz

    texturearray = CuTextureArray(p)
    # Create a texture object and bind it to the texture memory created above
    # p_texture = CuTexture(texturearray)
    p_texture = CuTexture(texturearray; filter_mode=CUDA.CU_TR_FILTER_MODE_LINEAR)

    # for iz=1:block:nangles
        @cuda threads=nthreads blocks=nblocks bp_parallel3d_kernel!(img, p_texture, size(p,1), size(p,2), nangles, nx, ny, nz, cons, scale)
    # end
end
