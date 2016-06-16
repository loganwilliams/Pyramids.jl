using Pyramids
using Base.Test, Images

function end_to_end(T, convert_to_arr=true)
    zone_plate = load("zone_plate.png")

    if convert_to_arr
        zone_plate = convert(Array, zone_plate)
        zone_plate = convert(Array, separate(zone_plate))
    end

    if typeof(T) <: ComplexSteerablePyramid
        scale = 0.5^(1/4);
        pyramid = ImagePyramid(zone_plate, T, scale=scale, num_orientations=8, max_levels=23, min_size=15, twidth=1.0)
    else
        pyramid = ImagePyramid(zone_plate, T, max_levels=23, min_size=15)
    end

    zone_plate_recon = toimage(pyramid)

    return all(zone_plate_recon .- zone_plate .< 1e-4)
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
