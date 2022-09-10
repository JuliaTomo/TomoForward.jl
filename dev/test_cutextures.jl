using Images, ColorTypes, FixedPointNumbers
using TestImages
using CUDA

# Get the input image. Use RGBA to have 4 channels since CUDA textures can have only 1, 2 or 4 channels.
# img_ = RGBA{N0f8}.(testimage("lighthouse"))
# img_ = zeros(Float32,10,20,40,1)
# img = reinterpret(NTuple{4,UInt8}, img_)
img = zeros(Float32,10,20,30)

# Create a texture memory object (CUDA array) and initilaize it with the input image content (from host).
texturearray = CuTextureArray(img)

# Create a texture object and bind it to the texture memory created above
texture = CuTexture(texturearray)

# Define an image warping kernel
function warp(dst, texture)
    tid = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    I = CartesianIndices(dst)
    @inbounds if tid <= length(I)
        i,j = Tuple(I[tid])
        u = Float32(i-1) / Float32(size(dst, 1)-1)
        v = Float32(j-1) / Float32(size(dst, 2)-1)
        x = u + 0.02f0 * CUDA.sin(30v)
        y = v + 0.03f0 * CUDA.sin(20u)
        dst[i,j] = texture[x,y, 0.2f0]
    end
    return
end

function configurator(kernel)
    config = launch_configuration(kernel.fun)

    threads = Base.min(length(outimg_d), config.threads)
    blocks = cld(length(outimg_d), threads)

    return (threads=threads, blocks=blocks)
end

# Create a 500x1000 CuArray for the output (warped) image
outimg_d = CuArray{eltype(img)}(undef, 500, 1000)

# Execute the kernel
@cuda config=configurator warp(view(outimg_d, 1:200,:), texture)
# @cuda threads = (size(outimg_d, 1), 1) blocks = (1, size(outimg_d, 2)) warp(outimg_d, texture)

# Get the output image into host memory and save it to a file
outimg = Array(outimg_d)
