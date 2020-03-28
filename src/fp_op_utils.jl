"forward project 3d reconstruction [ny nx nz] slice by slice onto phat [nangles detcount nz]"
function fp_2d_slices!(phat, A, u)
    @assert size(phat, 2) == size(u, 3)
    for slice=1:size(u, 3)
        bb = vec(view(phat, :, slice, :))
        phat[:,slice,:] = reshape( A*vec(u[:,:,slice]), size(phat, 1), size(phat, 3))
    end
end
