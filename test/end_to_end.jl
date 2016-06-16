function end_to_end()
    zone_plate = load("test/zone_plate.png")
    zone_plate = convert(Array, zone_plate)
    zone_plate = convert(Array, separate(zone_plate))

    scale = 0.5^(1/4);
    pyramid = ImagePyramid(zone_plate, scale, num_orientations=8, max_levels=23, min_size=15, twidth=1.0)

    zone_plate_recon = toimage(pyramid)

    return all(zone_plate_recon .- zone_plate .< 1e-4)
end