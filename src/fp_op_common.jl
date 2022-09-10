module fp_op_common
export SP, add_weight, check_vol_geom

struct SP{T<:AbstractFloat}
    I::Array{Int,1}
    J::Array{Int,1}
    V::Array{T,1}
    nrows::Int
    ncols::Int
end

# "Cons"
function SP(nrows, ncols)
    I = Array{Int,1}()
    J = Array{Int,1}()
    V = Array{Float64,1}()

    sz = Int(floor(nrows/100))*ncols
    sizehint!(I, sz)
    sizehint!(J, sz)
    sizehint!(V, sz)

    return SP(I, J, V, nrows, ncols)
end


"""
    add_weight()
    
Add the weight for sparse matrix A 
"""
function add_weight(A, iray, ivol, weight)
    push!(A.I, iray)
    push!(A.J, ivol)
    push!(A.V, weight)
    
    
    if iray > A.nrows || ivol > A.ncols
        println(iray, " ", ivol, " ", weight)
        error("index errors")
    end
end

function check_vol_geom(proj_geom, vol_width )
    if proj_geom.DetectorSpacingX > 0.0
        detwidth = (proj_geom.DetectorColCount*proj_geom.DetectorSpacingX)
        if vol_width > detwidth*10
            error("! Volume is too large compared to projection geometry. Please specify minX, maxX, minY, maxY when you call fp_op_parallel2d_strip")
        end
    end
end

end
