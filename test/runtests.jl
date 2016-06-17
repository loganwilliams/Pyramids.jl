using Pyramids
using Base.Test, Images

function end_to_end(T; convert_to_arr=true)
    test_im = rand(64, 64)

    if !convert_to_arr
        test_im = convert(Image, test_im)
    end

    if typeof(T) <: ComplexSteerablePyramid
        scale = 0.5^(1/4);
        pyramid = ImagePyramid(test_im, T, scale=scale, num_orientations=8, max_levels=23, min_size=15, twidth=1.0)
    else
        pyramid = ImagePyramid(test_im, T, max_levels=23, min_size=15)
    end

    test_im_recon = toimage(pyramid)

    return all((test_im_recon .- test_im) .< 0.0002)
end

println("Running end-to-end image comparison test for Complex Steerable Pyramid.")
@test end_to_end(ComplexSteerablePyramid())
println("\twithout array conversion")
@test end_to_end(ComplexSteerablePyramid(), convert_to_arr=false)

println("Running end-to-end image comparison test for Gaussian Pyramid.")
@test end_to_end(GaussianPyramid())
println("\twithout array conversion")
@test end_to_end(GaussianPyramid(), convert_to_arr=false)

println("Running end-to-end image comparison test for Laplacian Pyramid.")
@test end_to_end(LaplacianPyramid())

println("Tests passing.")
