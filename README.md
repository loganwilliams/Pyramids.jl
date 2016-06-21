# Pyramids.jl

A multi-scale image pyramid representation library for Julia.

Currently, Pyramids.jl can create, manipulate, and reduce Gaussian pyramids, Laplacian pyramids, and complex steerable pyramids. These image representations are used for a wide variety of computer vision and computational photography algorithms.

## Project Status

[![Build Status](https://travis-ci.org/loganwilliams/Pyramids.jl.svg?branch=master)](https://travis-ci.org/loganwilliams/Pyramids.jl) [![Coverage Status](https://coveralls.io/repos/github/loganwilliams/Pyramids.jl/badge.svg?branch=master)](https://coveralls.io/github/loganwilliams/Pyramids.jl?branch=master)

## Documentation

### Installation

Pyramids is not yet a registered Julia package. To install it, run the following in Julia:

    Pkg.clone("https://github.com/loganwilliams/Pyramids.jl.git")
    using Pyramids

### Usage

While parts of the Pyramids library are adapted from MatlabPyrTools, use of the library is quite different.

The Pyramids library is used through the `ImagePyramid` class. To create a new `ImagePyramid`, there a several possible constructors. Most typically, they take the form `ImagePyramid(im::Image, t::PyramidType)`.

Subtypes of `PyramidType` include:
 * The abstract type `SimplePyramid`
     - `GaussianPyramid`
     - `LaplacianPyramid`
 * `ComplexSteerablePyramid`

Additional parameters include the pyramid scale/slope, the maximum number of levels, number of orientation bands, and the minimum level size, as applicable. For a complete listing of parameters, view the source code.

Below is an example of loading an image and converting it to an `ImagePyramid`.

    using Images, Pyramids

    im = load("test_im.png")
    pyramid = ImagePyramid(im, ComplexSteerablePyramid(), scale=0.75, max_levels=23)

An `ImagePyramid` may also be created directly from an `Array`.

To manipulate an `ImagePyramid`, the `subband` and `update_subband` functions may be used.

    hi_residual = subband(pyramid, 0)
    pyramid_inverted_hi = update_subband(pyramid, 0, -hi_residual)

`update_subband` also has a mutating variant, `update_subband!`.

To convert an `ImagePyramid` back to its original image representation, the `toimage` function may be used.

    image_inverted_hi = toimage(pyramid_inverted_hi)

### Further examples

The `examples/` directory provides a further example of using Pyramids to implement "Phase-Based Frame Interpolation for Video," by Meyer et. al., from CVPR 2015. [1]

[1] https://www.disneyresearch.com/publication/phasebased/
