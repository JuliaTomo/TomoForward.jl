using Images, TestImages, ColorTypes, FixedPointNumbers
using CuArrays, CUDAnative
using CuTextures

# Get the input image. Use RGBA to have 4 channels since CUDA textures can have only 1, 2 or 4 channels.
#img = RGBA{N0f8}.(testimage("lighthouse"))

img = zeros(Float32,10,20,20)

# Create a texture memory object (CUDA array) and initilaize it with the input image content (from host).
texturearray = CuTextureArray(img)

# Create a texture object and bind it to the texture memory created above
texture = CuTexture(texturearray)

# Define an image warping kernel
function warp(dst, texture)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    u = (Float32(i) - 1f0) / (Float32(size(dst, 1)) - 1f0)
    v = (Float32(j) - 1f0) / (Float32(size(dst, 2)) - 1f0)
    x = u + 0.02f0 * CUDAnative.sin(30v)
    y = v + 0.03f0 * CUDAnative.sin(20u)
    @inbounds dst[i,j] = texture(x,y)
    return
end

# Create a 500x1000 CuArray for the output (warped) image
outimg_d = CuArray{eltype(img)}(undef, 500, 1000)

# Execute the kernel
@cuda threads = (size(outimg_d, 1), 1) blocks = (1, size(outimg_d, 2)) warp(outimg_d, texture)

# Get the output image into host memory and save it to a file
outimg = Array(outimg_d)
