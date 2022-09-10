using SparseArrays

using .fp_op_common

function fp_op_parallel2d_strip(proj_geom, H, W, mask_exclude=nothing)
    minX = -W // 2
    maxX = +W // 2
    minY = -H // 2
    maxY = +H // 2
    
    check_vol_geom(proj_geom, maxX-minX )
    fp_op_parallel2d_strip(proj_geom, H, W, minX, maxX, minY, maxY, mask_exclude)
end

function fp_op_parallel2d_strip(proj_geom::ProjGeom, vol_geom::VolGeom, mask_exclude=nothing)
    fp_op_parallel2d_strip(proj_geom, vol_geom.ny, vol_geom.nx, vol_geom.minX, vol_geom.maxX, vol_geom.minY, vol_geom.maxY, mask_exclude)
end


# The following code was ported to Julia from
# https://github.com/astra-toolbox/astra-toolbox/blob/master/include/astra/ParallelBeamStripKernelProjector2D.inl
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
function fp_op_parallel2d_strip(proj_geom::ProjGeom, H, W, minX, maxX, minY, maxY, mask_exclude=nothing)

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
    
    pixel_area = pixelspacingX * pixelspacingY

    # for each angle
    for i in 1:nangles
        vector = proj_geom.Vectors[i,:]

        # ray vector (D+aR)
        rx = vector[1]
        ry = vector[2]
        
        vertical = abs(rx) < abs(ry)

        Dx0 = -rx * src_origin - vector[5] * detcount / 2.
        Dy0 = -ry * src_origin - vector[6] * detcount / 2.

        raywidth = abs(vector[5]*ry - vector[6]*rx) / sqrt(rx*rx+ry*ry)
        rel_pixel_area = pixel_area / raywidth
        
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
            
            DLx = Dx0 + (j-0.5) * vector[5]
            DLy = Dy0 + (j-0.5) * vector[6]
            DRx = DLx + vector[5]
            DRy = DLy + vector[6]
            
            isin = false
            
            if (vertical)
                rxry = rx / ry
                Δc = - pixelspacingY * rxry * invpixelspacingX

                S = 0.5 - 0.5 * abs(rxry)
                T = 0.5 + 0.5 * abs(rxry)

                invTS = 1.0 / (T - S)

                # D+aR and E+cF
                cL = (DLx + (Ey - DLy)*rxry - Ex) * invpixelspacingX
                cR = (DRx + (Ey - DRy)*rxry - Ex) * invpixelspacingX
                
                if cR < cL
                    tmp_ = cL
                    cL = cR
                    cR = tmp_
                end

                # for each row
                for row = 0:H-1
                    colL = Int(floor(cL-0.5+S))
                    colR = Int(floor(cR+1.5-S))
                    colL = max(0, colL)
                    colR = min(colR, W-1)
                    
                    tmp = colL
                    offsetL = cL - tmp
                    offsetR = cR - tmp

                    for col=colL:colR
                        ivol = (col) * H + row+1
                        
                        weight = 0.0
                        
                        # right
                        if T <= offsetR
                            weight = 1.0
                        elseif S < offsetR
                            weight = 1.0 - 0.5*(T-offsetR)*(T-offsetR)*invTS
                        elseif -S < offsetR
                            weight = 0.5 + offsetR
                        elseif -T < offsetR
                            weight = 0.5*(offsetR+T)*(offsetR+T)*invTS
                        end
                        
                        # left
                        if T <= offsetL
                            weight -= 1.0
                        elseif S < offsetL
                            weight -= 1.0 - 0.5*(T-offsetL)*(T-offsetL)*invTS
                        elseif -S < offsetL
                            weight -= 0.5 + offsetL
                        elseif -T < offsetL
                            weight -= 0.5*(offsetL+T)*(offsetL+T)*invTS
                        end
                        
                        # if on_the_fly == false
                            add_weight(A, iray, ivol, weight*rel_pixel_area)
                        # else
                            # p[iray] += img[ivol] * weight
                        # end
                        offsetL -= 1.0
                        offsetR -= 1.0
                    end
                    cL += Δc
                    cR += Δc
                end 
            else
                ryrx = ry/rx
                Δr  = -pixelspacingX * ryrx * invpixelspacingY

                S = 0.5 - 0.5*abs(ryrx)
                T = 0.5 + 0.5*abs(ryrx)
                invTS = 1.0 / (T - S)

                rL = -(DLy + (Ex-DLx)*ryrx - Ey) * invpixelspacingY
                rR = -(DRy + (Ex-DRx)*ryrx - Ey) * invpixelspacingY
                
                if rR < rL
                    tmp_ = rL
                    rL = rR
                    rR = tmp_
                end

                for col=0:W-1
                    row_top = Int(floor(rL-0.5+S))
                    row_bot = Int(floor(rR+1.5-S))
                    
                    row_top = max(row_top, 0)
                    row_bot = min(row_bot, H-1)
                    
                    tmp = float(row_top)
                    offsetL = rL - tmp
                    offsetR = rR - tmp
                    
                    for row=row_top:row_bot
                        ivol = col * H + row+1
                        
                        weight = 0.0
                        
                        if T <= offsetR
                            weight = 1.0
                        elseif S < offsetR
                            weight =1.0 - 0.5*(T-offsetR)^2*invTS
                        elseif -S < offsetR
                            weight = 0.5 + offsetR
                        elseif -T < offsetR
                            weight = 0.5*(offsetR+T)*(offsetR+T)*invTS
                        end
                        
                        # left edge
                        if T <= offsetL
                            weight -= 1.0
                        elseif S < offsetL
                            weight -= 1.0-0.5*(T-offsetL)^2*invTS
                        elseif -S < offsetL
                            weight -= 0.5 + offsetL
                        elseif -T < offsetL
                            weight -= 0.5*(offsetL+T)^2*invTS
                        end
                        
                        # if on_the_fly == false
                            add_weight(A, iray, ivol, weight*rel_pixel_area)
                        # else
                            # p[iray] += img[ivol] * weight
                        # end
                        
                        offsetL -= 1.0
                        offsetR -= 1.0
                    end
                    
                    rL += Δr
                    rR += Δr
                end
            end
        end
    end
    
    A = sparse(A.I, A.J, A.V, A.nrows, A.ncols)
    dropzeros!(A)
    return A
end
