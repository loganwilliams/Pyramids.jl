# This file shows an example use of Pyramids.jl to implement an algorithm that requires Complex Steerable Pyramids.
#
# This is an implementation of "Phase-Based Frame Interpolation for Video," by Meyer et. al., from CVPR 2015. [1]
#
# [1] https://www.disneyresearch.com/publication/phasebased/

include("Pyramids.jl")

using Images
using ColorTypes
using Pyramids
using Interpolations

type PhasePyramid <: PyramidType
end

function shift_correction(phi_delta::ImagePyramid; shift_limit=0.5)
    corrected_phase_delta = ImagePyramid(phi_delta)

    for level = corrected_phase_delta.num_levels:-1:1
        for orientation = 1:corrected_phase_delta.num_orientations
            phi_l = subband(corrected_phase_delta, level, orientation=orientation)
            dims_l = size(phi_l)
            phi_limit = (π / (phi_delta.scale^(phi_delta.num_levels - level))) * shift_limit

            if level < phi_delta.num_levels
                phi_lp1 = subband(corrected_phase_delta, level+1, orientation=orientation)
                dims_lp1 = size(phi_lp1)
                phi_lp1_itp = interpolate(phi_lp1, BSpline(Linear()), OnGrid())
               
                for r = 1:dims_l[1]
                    for c = 1:dims_l[2]

                        r_lp1 = ((r-1)/(dims_l[1]-1)) * (dims_lp1[1]-1) + 1
                        c_lp1 = ((c-1)/(dims_l[2]-1)) * (dims_lp1[2]-1) + 1

                        # unwrap
                        while phi_l[r,c] > (phi_lp1_itp[r_lp1,c_lp1]/phi_delta.scale + π)
                            phi_l[r,c] -= 2*π
                        end

                        while phi_l[r,c] < (phi_lp1_itp[r_lp1,c_lp1]/phi_delta.scale - π)
                            phi_l[r,c] += 2*π
                        end

                        confidence = atan2(sin(phi_l[r,c] - phi_lp1_itp[r_lp1,c_lp1]/phi_delta.scale), cos(phi_l[r,c] - phi_lp1_itp[r_lp1,c_lp1]/phi_delta.scale))

                        if (abs(confidence) > π/2)
                            phi_l[r,c] = phi_lp1_itp[r_lp1,c_lp1]/phi_delta.scale
                        end
                        
                        if abs(phi_l[r,c]) > phi_limit
                            phi_l[r,c] = phi_lp1_itp[r_lp1,c_lp1]/phi_delta.scale
                        end

                    end
                end
            else
                phi_l[abs(phi_l) .> phi_limit] = 0
            end

            update_subband!(corrected_phase_delta, level, phi_l, orientation=orientation)
        end
    end

    return corrected_phase_delta
end

function adjust_phase(phase_delta::ImagePyramid, corrected_phase_delta::ImagePyramid)
    adjusted_phase = adjust_phase(phase_delta.pyr, corrected_phase_delta.pyr)
    return ImagePyramid(adjusted_phase, phase_delta.pind, phase_delta.scale, PhasePyramid())
end

function adjust_phase(phase_delta::Array, corrected_phase_delta)
    adjusted_phase_delta = copy(phase_delta)

    for n = 1:length(adjusted_phase_delta)
        while adjusted_phase_delta[n] > (corrected_phase_delta[n] + π)
            adjusted_phase_delta[n] -= 2*π
        end

        while adjusted_phase_delta[n] < (corrected_phase_delta[n] - π)
            adjusted_phase_delta[n] += 2*π
        end
    end

    return adjusted_phase_delta
end

# function joint_unwrap!(phase1, phase2)


# end

function blend_and_interpolate(pyramid1::ImagePyramid, pyramid2::ImagePyramid, phase_delta::ImagePyramid, alpha)
    A1 = abs(pyramid1.pyr)
    A2 = abs(pyramid2.pyr)
    phase1 = angle(pyramid1.pyr)
    phi_delta = phase_delta.pyr

    new_A = A1 * (1 - alpha) + A2 * alpha
    new_phase = phase1 + alpha * phi_delta
    new_pyr = new_A .* exp(complex(0,1) * new_phase)

    new_pyramid = ImagePyramid(new_pyr, pyramid1.pind, pyramid1.scale, ComplexSteerablePyramid())

    update_subband!(new_pyramid, 0, subband(pyramid2, 0)*0)

    if (alpha > 0.5)
        update_subband!(new_pyramid, 0, subband(pyramid2, 0))
    else
        update_subband!(new_pyramid, 0, subband(pyramid1, 0))
    end

    return new_pyramid

end

println("Loading images")

im1 = load("frame_0.jpg")
im1 = convert(Image{Lab}, float32(im1))

im2 = load("frame_1.jpg")
im2 = convert(Image{Lab}, float32(im2))

im1 = convert(Array, separate(im1))
im2 = convert(Array, separate(im2))

L1 = im1[:,:,1]
L2 = im2[:,:,1]

println("Converting images to complex steerable pyramids")

pyramid1 = ImagePyramid(L1, ComplexSteerablePyramid(), scale=0.5^0.25)
pyramid2 = ImagePyramid(L2, ComplexSteerablePyramid(), scale=0.5^0.25)

phase_delta = ImagePyramid(mod(angle(pyramid2.pyr) - angle(pyramid1.pyr) + π, 2*π) - π, pyramid1.pind, 0.5^0.25, PhasePyramid())
corrected_phase_delta = shift_correction(phase_delta, shift_limit=1)
adjusted_phase_delta = adjust_phase(phase_delta, corrected_phase_delta)

for alpha = 0:0.2:1.0
    println("Generating image with alpha $(alpha)")
    newpyr = blend_and_interpolate(pyramid1, pyramid2, adjusted_phase_delta, alpha)
    newim = toimage(newpyr)

    newLabIm = zeros(im1)
    newLabIm[:,:,1] = newim
    newLabIm[:,:,2] = (1 - alpha) * im1[:,:,2] + alpha * im2[:,:,2]
    newLabIm[:,:,3] = (1 - alpha) * im1[:,:,3] + alpha * im2[:,:,3]

    newLabIm = convert(Image, newLabIm)
    newLabIm.properties["colorspace"] = "Lab"
    newLabIm = convert(Image{RGB}, newLabIm)
    save("interpolated_frame_$(alpha).png", newLabIm')
end