using SparseArrays

using .fp_op_common

function fp_op_parallel2d_line(proj_geom::ProjGeom, H, W, mask_exclude=nothing)
    minX = -W // 2
    maxX = +W // 2
    minY = -H // 2
    maxY = +H // 2
    
    check_vol_geom(proj_geom, maxX-minX )
    fp_op_parallel2d_line(proj_geom, H, W, minX, maxX, minY, maxY, mask_exclude)
end

function fp_op_parallel2d_line(proj_geom::ProjGeom, vol_geom::VolGeom, mask_exclude=nothing)
    fp_op_parallel2d_line(proj_geom, vol_geom.ny, vol_geom.nx, vol_geom.minX, vol_geom.maxX, vol_geom.minY, vol_geom.maxY, mask_exclude)
end


# The following code was ported to Julia from
# https://github.com/astra-toolbox/astra-toolbox/blob/master/include/astra/ParallelBeamLinearKernelProjector2D.inl
# /*
# -----------------------------------------------------------------------
# Copyright: 2010-2018, imec Vision Lab, University of Antwerp
#            2014-2018, CWI, Amsterdam
# Contact: astra@astra-toolbox.com
# Website: http://www.astra-toolbox.com/
# This file is part of the ASTRA Toolbox.
# The ASTRA Toolbox is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# The ASTRA Toolbox is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with the ASTRA Toolbox. If not, see <http://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------
# */

"""
    get_sys_parallel2d

Compute the sparse system matrix
for the image size of [H x W]

# optional arguments
img ::Array{AbstractFloat,1} vectorized image
"""
function fp_op_parallel2d_line(proj_geom::ProjGeom, H, W, minX, maxX, minY, maxY, mask_exclude=nothing)
    # Partially ported from ASTRA-Toolbox
    # https://github.com/astra-toolbox/astra-toolbox/blob/master/include/astra/ParallelBeamLineKernelProjector2D.inl

    nangles = size(proj_geom.Vectors, 1)
    detcount = proj_geom.DetectorColCount
    
    pixelspacingX = (maxX-minX) / W
    pixelspacingY = (maxY-minY) / H
    
    invpixelspacingX = 1. / pixelspacingX
    invpixelspacingY = 1. / pixelspacingY

    src_origin = 0.0

    Ex = minX + pixelspacingX*0.5 # min x
    Ey = maxY - pixelspacingY*0.5 # max y

    
    A = SP(nangles*detcount, H*W)
    
    # for each angle
    for i in 1:nangles
        vector = proj_geom.Vectors[i,:]

        # ray vector (D+aR)
        rx = vector[1]
        ry = vector[2]
        
        vertical = abs(rx) < abs(ry)

        Dx0 = -rx * src_origin - vector[5] * detcount / 2.
        Dy0 = -ry * src_origin - vector[6] * detcount / 2.
        
        # for each ray
        for j = 1:detcount
            if ~isnothing(mask_exclude)
                # (for sinogram inpainting)
                # if a mask_exclude is given, check if the pixel is not considered
                if mask_exclude[i,j] == 1
                    continue
                end
            end

            iray = (j-1)*nangles + i

            Dx = Dx0 + (j-0.5) * vector[5]
            Dy = Dy0 + (j-0.5) * vector[6]

            # E: left top position of an image
            # R: ray direction vector
            # D: detector position
            
            isin = false

            if (vertical)
                rxry = rx / ry
                len = pixelspacingX * sqrt( rx*rx + ry*ry ) / abs(ry)
                Δc = - pixelspacingY * rxry * invpixelspacingX

                S = 0.5 - 0.5 * abs(rxry)
                T = 0.5 + 0.5 * abs(rxry)

                invLTS = len / (T - S)

                # D+aR and E+cF
                c = (Dx + (Ey - Dy)*rxry - Ex) * invpixelspacingX

                # for each row
                for row = 0:H-1
                    col = Int(floor(c+0.5))
                    if (col < -1 || col > W)
                        if (isin == false)
                            c += Δc
                            continue
                        else break; end
                    end

                    offset = c - float(col)

                    # left
                    if (offset < -S)
                        weight = (offset + T) * invLTS
                        # zero-based: col-1, row
                        ivol = (col-1) * H + (row+1)
                        if col > 0
                                add_weight(A, iray, ivol, len-weight)
                        end

                        ivol = (col+1-1) * H + (row+1)
                        if col >= 0 && col < W
                                add_weight(A, iray, ivol, weight)
                        end

                    elseif (S < offset)
                        weight = (offset - S) * invLTS

                        ivol = (col) * H + (row+1)
                        if col >= 0 && col < W
                                add_weight(A, iray, ivol, len-weight)
                        end

                        ivol = (col+1) * H + (row+1)
                        if (col+1 < W)
                                add_weight(A, iray, ivol, weight)
                        end

                    elseif (col >= 0 && col < W)
                        ivol = col * H + (row+1)
                            add_weight(A, iray, ivol, len)
                    end

                    c += Δc
                    isin = true
                end
            else
                ryrx = ry/rx
                len = pixelspacingY * sqrt(rx*rx+ry*ry) / abs(rx)
#                 len = pixelspacingX*pixelspacingY / raywidth
                Δr  = -pixelspacingX * ryrx * invpixelspacingY

                S = 0.5 - 0.5*abs(ryrx)
                T = 0.5 + 0.5*abs(ryrx)
                invLTS = len / (T - S)

                r = -(Dy + (Ex-Dx)*ryrx - Ey) * invpixelspacingY

                for col=0:W-1
                    
                    row = Int(floor(r+0.5))
                    if (row < -1 || row > H)
                        if isin == false
                            r += Δr
                            continue
                        else
                            break
                        end
                    end
                    
                    offset = r - float(row)

                    if (offset < -S)
                        weight = (offset + T) * invLTS
                        
                        ivol = col * H + (row-1) + 1   
                        if row > 0
                            add_weight(A, iray, ivol, len-weight)
                        end

                        ivol = col * H + (row) + 1   
                        if (row >= 0 && row < H)
                            add_weight(A, iray, ivol, weight)
                        end

                    elseif (S < offset)
                        weight = (offset - S) * invLTS

                        ivol = col * H + (row+1)
                        if row >= 0 && row < H
                            add_weight(A, iray, ivol, len-weight)
                        end

                        # zero-based: col, row+1
                        ivol = col * H + (row+1) + 1
                        
                        if row+1 < H
                            add_weight(A, iray, ivol, weight)
                        end

                    elseif row >= 0 && row < H
                        ivol = col * H + (row+1)
                        add_weight(A, iray, ivol, len)
                    end

                    isin = true
                    r += Δr
                end
            end
        end
    end
    
    A = sparse(A.I, A.J, A.V, A.nrows, A.ncols)
    dropzeros!(A)
    return A
end
