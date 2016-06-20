module Pyramids

using Images, Interpolations, Colors

export ImagePyramid, PyramidType, ComplexSteerablePyramid, LaplacianPyramid, GaussianPyramid
export subband, toimage, update_subband, update_subband!, test

"Abstract supertype for the variety of pyramids that can be constructed."
abstract PyramidType
abstract SimplePyramid <: PyramidType

"""Type that indicates a complex steerable pyramid. [1]

[1] http://www.cns.nyu.edu/~eero/steerpyr/"""
type ComplexSteerablePyramid <: PyramidType
end

"""Type that indicates a Laplacian pyramid. [1]

[1] persci.mit.edu/pub_pdfs/pyramid83.pdf"""
type LaplacianPyramid <: SimplePyramid
end

"""Type that indicates a Gaussian pyramid. [1]

[1] http://persci.mit.edu/pub_pdfs/RCA84.pdf"""
type GaussianPyramid <: SimplePyramid
end

"""Type that represents a concrete pyramidal representation of a given input image. Each type of pyramid has its own parameters. The basic construction method is, for example

```pyramid = ImagePyramid(im, ComplexSteerablePyramid(), scale=0.5^0.25)```

See the code for more information on optional arguments."""
type ImagePyramid
    scale::Real
    pyr::Array
    pind::Array{Int}
    num_orientations::Int
    num_levels::Int
    t::PyramidType
    pyr_dict::Dict

    function ImagePyramid(pyr::ImagePyramid)
        return deepcopy(pyr)
    end

    function ImagePyramid(im::Array, t::ComplexSteerablePyramid; scale=0.5, min_size=15, num_orientations=8, max_levels=23, twidth=1)
        this = new()

        this.t = ComplexSteerablePyramid()
        this.scale = scale
        this.num_orientations = num_orientations

        im_dims = size(im)
        h = im_dims[1]
        w = im_dims[2]

        this.num_levels = min(ceil(log2(minimum([h w]))/log2(1/scale) - (log2(min_size)/log2(1/scale))),max_levels);

        pyr_dict, mtx, harmonics = build_complex_steerable_pyramid(im, this.num_levels, this.num_levels, order=num_orientations-1, twidth=twidth, scale=scale)

        # print(pyr_dict)

        # this.pyr = pyr
        # this.pind = pind
        this.pyr_dict = pyr_dict

        return this
    end

    function ImagePyramid(im::Image, t::ComplexSteerablePyramid; scale=0.5, min_size=15, num_orientations=8, max_levels=23, twidth=1)
        im_arr = image_to_array(im)

        return ImagePyramid(im_arr, t, scale=scale, min_size=min_size, num_orientations=num_orientations, max_levels=23, twidth=1)
    end

    function ImagePyramid(im::Array, t::GaussianPyramid; min_size=15, max_levels=23, filter=[0.0625; 0.25; 0.375; 0.25; 0.0625])
        this = new()

        pyr, pind, num_levels = generate_gaussian_pyramid(im, min_size=min_size, max_levels=max_levels, filter=filter)

        this.scale = 0.5
        this.pyr = pyr
        this.pind = pind
        this.num_orientations = 1
        this.num_levels = num_levels
        this.t = t

        return this
    end

    function ImagePyramid(im::Array, t::LaplacianPyramid; min_size=15, max_levels=23, filter=[0.0625; 0.25; 0.375; 0.25; 0.0625])
        this = new()

        pyr, pind, num_levels = generate_laplacian_pyramid(im, min_size=min_size, max_levels=max_levels, filter=filter)

        this.scale = 0.5
        this.pyr = pyr
        this.pind = pind
        this.num_orientations = 1
        this.num_levels = num_levels
        this.t = t

        return this
    end

    function ImagePyramid(im::Image, t::SimplePyramid; min_size=15, max_levels=23, filter=[0.0625; 0.25; 0.375; 0.25; 0.0625])
        im_arr = image_to_array(im)

        return ImagePyramid(im_arr, t, min_size=min_size, max_levels=max_levels, filter=filter)
    end

    function ImagePyramid(pyr::Array, pind::Array{Int}, scale, t)
        this = new()

        this.scale = scale
        this.pyr = pyr
        this.pind = pind
        this.num_orientations = num_orientation_bands(pind)
        this.num_levels = pyramid_height(pind)
        this.t = t

        return this
    end

end

##############################
# Functions

"""
    subband(pyramid::ImagePyramid, level; orientation=())

Returns the sub-band of the `pyramid` at `level` and at `orientation.` If orientation is not provided, `subband` assumes that the pyramid is not oriented.

Level 0 is always the high frequency residual."""
function subband(pyramid::ImagePyramid, level; orientation = ())
    if isempty(orientation)
        band = level + 1
    else
        band = (level-1)*pyramid.num_orientations + orientation + 2
    end

    return pyramid_subband(pyramid.pyr, pyramid.pind, band)
end

function update_subband(pyramid::ImagePyramid, level, new_subband; orientation = ())
    newpyramid = ImagePyramid(pyramid)

    if isempty(orientation)
        band = level + 1
    else
        band = (level-1)*newpyramid.num_orientations + orientation + 2
    end

    newpyr = update_pyramid_subband(newpyramid.pyr, newpyramid.pind, band, new_subband)
    newpyramid.pyr = newpyr

    return newpyramid
end

function update_subband!(pyramid::ImagePyramid, level, new_subband; orientation = ())
    if isempty(orientation)
        band = level + 1
    else
        band = (level-1)*pyramid.num_orientations + orientation + 2
    end

    pyramid.pyr = update_pyramid_subband(pyramid.pyr, pyramid.pind, band, new_subband)
    
    return pyramid
end

"Converts a pyramid to a 2-D array (not of type `Image`!)"
function toimage(pyramid::ImagePyramid)
    if typeof(pyramid.t) <: ComplexSteerablePyramid
        im = reconstruct_complex_steerable_pyramid(pyramid.pyr, pyramid.pind, scale=pyramid.scale)
    elseif typeof(pyramid.t) <: LaplacianPyramid
        im = reconstruct_laplacian_pyramid(pyramid)
    elseif typeof(pyramid.t) <: GaussianPyramid
        im = subband(pyramid, 0)
    else
        error("Unsupported pyramid type $(typeof(pyramid.T))")
    end

    return im
end

##############################
# Private helper functions.

function image_to_array(im::Image)
    # Please someone tell me a better way of checking for image type.
    if typeof(im.data[1]) <: Gray
        im_arr = separate(convert(Array{Float64,2}, im))
    else
        warn("ImagePyramid only supports 2D images currently. Returning luminance channel.")
        im_lab = convert(Image{Lab}, im)
        im_arr = convert(Array{Float64,3}, separate(im))[:,:,1] 
    end

    return im_arr
end

##############################
# Private functions for building Gaussian and Laplacian pyramids.

function convolve_reflect(im, filter)
    filter_len = length(filter)
    filter_offset::Int = (filter_len-1)/2

    padded_im = zeros((size(im,1) + filter_len*2), (size(im,1) + filter_len*2))

    padded_im[(filter_len+1):(end-filter_len), (filter_len+1):(end-filter_len)] = im
    
    padded_im[(filter_len+1):(end-filter_len), 1:filter_len] = flipdim(padded_im[(filter_len+1):(end-filter_len), (filter_len+2):(2*filter_len+1)], 2)
    padded_im[(filter_len+1):(end-filter_len), (end-filter_len+1):end] = flipdim(padded_im[(filter_len+1):(end-filter_len), (end-2*filter_len):(end-filter_len-1)], 2)
    padded_im[1:filter_len, (filter_len+1):(end-filter_len),] = flipdim(padded_im[(filter_len+2):(1+2*filter_len), (filter_len+1):(end-filter_len)], 1)
    padded_im[(end-filter_len+1):end, (filter_len+1):(end-filter_len)] = flipdim(padded_im[(end-2*filter_len):(end-filter_len-1), (filter_len+1):(end-filter_len)], 1)

    new_im = conv2(filter, filter, padded_im)

    new_im = new_im[(1+filter_len+filter_offset):(end-filter_len-filter_offset), (1+filter_len+filter_offset):(end-filter_len-filter_offset)]
    return new_im
end

function reduce(im; filter=[0.0625; 0.25; 0.375; 0.25; 0.0625])
    new_im = convolve_reflect(im, filter)
    return new_im[1:2:end, 1:2:end]
end

# TODO: fix the expand to upsample and filter
function expand(im; filter=[0.0625; 0.25; 0.375; 0.25; 0.0625])
    new_im = zeros((size(im, 1)*2, size(im, 2)*2))
    new_im[1:2:end, 1:2:end] = im

    new_im = convolve_reflect(new_im, filter)
    return new_im
end

function generate_gaussian_pyramid(im; min_size=15, max_levels=23, filter=[0.0625; 0.25; 0.375; 0.25; 0.0625])
    im = convert(Array{Float64}, copy(im))

    pind = collect(size(im))';
    pyr = [];

    im_dims = collect(size(im))

    num_levels = min(max_levels, ceil(Int, log2(minimum(im_dims)) - log2(min_size)))

    for i = 1:num_levels
        # this is a weird hack to make the pyramid index in the same format as MatlabPyrTools.
        # TODO: change index format to something better suited for Julia (dictionary?)

        if i != 1
            pind = [pind; collect(size(im))'];
        end

        pyr = [pyr; im[:]]
        im = reduce(im, filter=filter)
    end

    return (pyr, pind, num_levels)
end

function generate_laplacian_pyramid(im; min_size=15, max_levels=23, filter=[0.0625; 0.25; 0.375; 0.25; 0.0625])
    im = convert(Array{Float64}, copy(im))

    pind = collect(size(im))';
    pyr = [];

    im_dims = collect(size(im))

    num_levels = min(max_levels, ceil(Int, log2(minimum(im_dims)) - log2(min_size)))

    filter_offset::Int = (length(filter)-1)/2

    for i = 1:num_levels
        reduced_im = reduce(im, filter=filter)
        next_im = expand(reduced_im)
        diff_im = im - next_im

        pyr = [pyr; diff_im[:]]

        if i != 1
            pind = [pind; collect(size(diff_im))']
        end

        im = reduced_im
    end

    pyr = [pyr; im[:]]
    pind = [pind; collect(size(im))']

    return (pyr, pind, num_levels)
end

function reconstruct_laplacian_pyramid(pyramid::ImagePyramid)
    output_im::Array{Float64} = subband(pyramid, pyramid.num_levels)

    for i = (pyramid.num_levels-1):-1:0
        output_im = expand(output_im)
        output_im += subband(pyramid, i)
    end

    return output_im
end

##############################
# Private functions for building complex steerable pyramids.
# Adapted from MatlabPyrTools. (https://github.com/LabForComputationalVision/matlabPyrTools)
# This all runs significantly slower than it should, particularly the pyramid reconstruction.

function construct_steering_matrix(harmonics, angles; even = true)
    numh = 2*length(harmonics) - any(harmonics .== 0)

    imtx = zeros(length(angles), numh)
    col = 1
    for h in harmonics
        args = h * angles

        if h == 0
            imtx[:,col] = ones(angles)
            col += 1
        elseif !even
            imtx[:,col] = sin(args)
            imtx[:,col+1] = -cos(args)
            col += 2
        else
            imtx[:,col] = cos(args)
            imtx[:,col+1] = sin(args)
            col += 2
        end
    end

    r = rank(imtx)

    if (r != numh) && (r != length(angles))
        warning("Matrix is not full rank")
    end

    return pinv(imtx)
end

function raisedcosine(width=1, position=0, values=[0,1])
    sz = 256 # arbitrary
    X = pi * (-sz-1:1) / (2 * sz)
    Y = values[1] + (values[2]-values[1]) * cos(X).^2

    Y[1] = Y[2]
    Y[sz+3] = Y[sz+2]
    X = position + (2*width/pi) * (X + pi/4)

    return (X, Y)
end

function build_complex_steerable_pyramid_level(lodft, log_rad, Xrcos, Yrcos, angleI, ht, nbands, scale, nScales, im_dims)    
    if ht <= 0
        
    else
        

        pyr, pind = build_complex_steerable_pyramid_level(lodft, log_rad, Xrcos, Yrcos, angleI, ht-1, nbands, scale, nScales, im_dims)

        
       
    end

    return (pyr, pind)
end

function build_complex_steerable_pyramid(im, height, nScales; order=3, twidth=1, scale=0.5)
    pyramid_bands = Dict{Integer, Union{Array, Dict{Integer, Array}}}()

    num_orientations = order + 1

    if mod(num_orientations, 2) == 0
        harmonics = ((0:((num_orientations/2)-1))*2 + 1)'
    else
        harmonics = ((0:((num_orientations-1)/2))*2)'
    end

    steeringmatrix = construct_steering_matrix(harmonics, pi*(0:(num_orientations-1))/num_orientations, even=true)

    im_dims = collect(size(im))
    ctr = ceil(Int, (im_dims+0.5)/2)

    angle = broadcast(atan2, (((1:im_dims[1]) - ctr[1]) ./ (im_dims[1]/2)), (((1:im_dims[2]) - ctr[2]) ./ (im_dims[2]/2))')
    log_rad = broadcast((x,y) -> log2(sqrt(x.^2 + y.^2)), (((1:im_dims[1]) - ctr[1]) ./ (im_dims[1]/2)), (((1:im_dims[2]) - ctr[2]) ./ (im_dims[2]/2))')
    log_rad[ctr[1], ctr[2]] = log_rad[ctr[1], ctr[2]-1]
    log_rad0 = copy(log_rad)

    Xrcos, Yrcos = raisedcosine(twidth, (-twidth/2), [0 1])
    Xrcos0 = Xrcos

    Yrcos = sqrt(Yrcos)

    YIrcos = sqrt(1 - Yrcos.^2)

    YIrcosinterpolant = interpolate((Xrcos,), YIrcos, Gridded(Linear()))
    lo0mask = reshape(YIrcosinterpolant[log_rad], size(log_rad))

    imdft = fftshift(fft(im))

    lo0dft = imdft .* lo0mask

    pyr = []
    pind = []

    for ht = height:-1:1
        Xrcos = Xrcos - log2(1/scale)
        lutsize = 1024
        Xcosn = pi*(-(2*lutsize+1):(lutsize+1))/lutsize
        order = num_orientations - 1

        cnst = (2^(2*order))*(factorial(order)^2)/(num_orientations*factorial(2*order))

        alfa = mod(pi+Xcosn, 2*pi) - pi
        Ycosn = 2*sqrt(cnst) * (cos(Xcosn).^order) .* (abs(alfa) .< pi/2)

        Yrcosinterpolant = interpolate((Xrcos,), Yrcos, Gridded(Linear()))
        himask = reshape(Yrcosinterpolant[log_rad], size(log_rad))

        pyramid_level = Dict{Integer, Array}()

        for b in 1:num_orientations
            Ycosninterpolant = interpolate((Xcosn + pi*(b-1)/num_orientations,), Ycosn, Gridded(Linear()))
            anglemask = reshape(Ycosninterpolant[angle], size(angle))
            banddft = ((complex(0,-1)).^(num_orientations-1)) .* lo0dft .* anglemask .* himask
            pyramid_level[b] = ifft(ifftshift(banddft))
        end

        pyramid_bands[height-ht+1] = pyramid_level

        dims = collect(size(lo0dft))
        ctr = ceil(Int, (dims+0.5)/2)

        lodims = round(Int, im_dims[1:2]*scale^(nScales-ht+1))

        loctr = ceil(Int, (lodims+0.5)/2)
        lostart = ctr - loctr+1
        loend = lostart + lodims - 1

        log_rad = log_rad[lostart[1]:loend[1], lostart[2]:loend[2]]
        angle = angle[lostart[1]:loend[1], lostart[2]:loend[2]]
        lodft = lo0dft[lostart[1]:loend[1], lostart[2]:loend[2]]
        YIrcos = abs(sqrt(1 - Yrcos.^2))

        YIrcosinterpolant = interpolate((Xrcos,), YIrcos, Gridded(Linear()))
        lomask = reshape(YIrcosinterpolant[log_rad], size(log_rad))

        lo0dft = lomask .* lodft
    end

    pyramid_bands[height+1] = real(ifft(ifftshift(lo0dft)))

    Yrcosinterpolant = interpolate((Xrcos0,), Yrcos, Gridded(Linear()))

    hi0mask = reshape(Yrcosinterpolant[log_rad0], size(log_rad0))
    hi0dft =  imdft .* hi0mask;
    pyramid_bands[0] = ifft(ifftshift(hi0dft));

    return (pyramid_bands, steeringmatrix, harmonics)
end

function pyramid_subband_index(pind, band)
    if (band > size(pind,1)) || (band < 1)
        error("band must be between 1 and number of pyramid bands ($(size(pind,1))).")
    end

    if size(pind, 2) != 2
        error("pind must be an Nx2 matrix indicating the size of the pyramid subbands")
    end

    ind = 1

    for l in 1:(band-1)
        ind = ind + prod(pind[l,:])
    end

    return round(Int, ind):round(Int, (ind + prod(pind[band,:]) - 1))
end

function pyramid_subband(pyr, pind, band)
    return reshape( copy(pyr[pyramid_subband_index(pind, band)]), (round(Int, pind[band,1]), round(Int, pind[band,2])))
end

function update_pyramid_subband(pyr, pind, band, newband)
    newpyr = copy(pyr)
    index = pyramid_subband_index(pind, band)
    newpyr[index] = newband
    return newpyr
end

function make_angle_grid(sz, phase=0, origin=-1)
    if length(sz) == 1
        sz = [sz, sz]
    end

    if origin == -1
        origin = (sz + 1)/2
    end

    xramp = ones(round(Int, sz[1]), 1) * collect((1:sz[2]) - origin[2])'
    yramp = collect((1:sz[1]) - origin[1]) * ones(1, round(Int, sz[2]))

    res = atan2(yramp, xramp)

    res = mod(res+(pi-phase), 2*pi) - pi

    return res
end

function num_orientation_bands(pind)
    if size(pind,1) == 2
        nbands = 0
    else
        b = 3
        while (b <= size(pind,1)) && all(pind[b,:] == pind[2,:])
            b += 1
        end

        nbands = b - 2
    end

    return nbands
end

function pyramid_height(pind)
    nbands = num_orientation_bands(pind)

    if size(pind,1) > 2
        ht = (size(pind,1) - 2) / nbands
    else
        ht = 0
    end

    return convert(Int, ht)
end

function reconstruct_steerable_pyramid_level(pyr, pind, log_rad, Xrcos, Yrcos, angle, nbands, levs, bands, scale)
    lo_ind = nbands + 1
    dims = convert(Array{Int}, pind[1,:])
    ctr = ceil(Int, (dims+0.5)/2)

    Xrcos = Xrcos - log2(1/scale)

    if any(levs .> 1)
        lodims = pind[1+bands[end],:]
        loctr = ceil(Int, (lodims+0.5)/2)
        lostart = ctr - loctr + 1
        loend = lostart + lodims - 1
        nlog_rad = log_rad[lostart[1]:loend[1], lostart[2]:loend[2]]
        nangle = angle[lostart[1]:loend[1], lostart[2]:loend[2]]
        
        if size(pind,1) > lo_ind
            subpyr = convert(Int, 1+sum(prod(pind[1:(lo_ind-1),:],2))):convert(Int, size(pyr, 1))
            nresdft = reconstruct_steerable_pyramid_level(pyr[subpyr], pind[lo_ind:size(pind,1),:], nlog_rad, Xrcos, Yrcos, nangle, nbands, levs-1, bands, scale)
        else
            nresdft = fftshift(fft(pyramid_subband(pyr, pind, lo_ind)))
        end

        YIrcos = sqrt(abs(1 - Yrcos.^2))
        YIrcosinterpolant = interpolate((Xrcos,), YIrcos, Gridded(Linear()))
        lomask = reshape(YIrcosinterpolant[nlog_rad], size(nlog_rad))

        resdft = complex(zeros((dims[1], dims[2])))
        resdft[lostart[1]:loend[1], lostart[2]:loend[2]] = nresdft .* lomask
    else
        resdft = complex(zeros((dims[1], dims[2])))
    end
    
    if any(levs .== 1)
        lutsize = 1024
        Xcosn = pi*(-(2*lutsize+1):(lutsize+1))/lutsize;
        order = nbands-1

        cnst = (2^(2*order))*(factorial(order)^2)/(nbands*factorial(2*order))
        Ycosn = sqrt(cnst) * (cos(Xcosn)).^order

        Yrcosinterpolant = interpolate((Xrcos,), Yrcos, Gridded(Linear()))
        himask = reshape(Yrcosinterpolant[log_rad], size(log_rad))

        ind = 1
        for b in 1:nbands
            if any(bands .== b)
                Ycosninterpolant = interpolate((Xcosn + pi*(b-1)/nbands,), Ycosn, Gridded(Linear()))
                anglemask = reshape(Ycosninterpolant[angle], size(angle))
                band = reshape(pyr[ind:ind+prod(dims)-1], (dims[1], dims[2]))
                banddft = fftshift(fft(band))
                resdft = resdft + complex(0,1)^(nbands-1) * banddft .* anglemask .* himask
            end
                
            ind = ind + prod(dims)
        end
    end

    return resdft
end

function sub_matrix(vec, sz, offset=1)
    vec = vec[:]
    sz = sz[:]

    if size(sz,1) != 2
        error("Dimensions must be a 2-vector")
    end

    return reshape(vec[offset:offset+prod(sz)-1], (convert(Int, sz[1]), (convert(Int, sz[2]))))
end


function reconstruct_steerable_pyramid(pyr, pind; levs="all", bands="all", twidth=1, scale=0.5)
    nbands = num_orientation_bands(pind)

    maxLev = 1 + pyramid_height(pind)

    if levs == "all"
        levs = collect(0:maxLev)'
    else
        if any(levs .> maxLev) || any(levs .< 0)
            error("Level numbers must be in the range [0, $(maxLev)]")
        end

        levs = collect(levs[:])
    end

    if bands == "all"
        bands = collect(1:nbands)'
    else
        if any(bands .< 1) || any(bands .> nbands)
            error("Band numbers must be in the range [1, $(nbands)]")
        end

        bands = collect(bands[:])
    end

    dims = pind[1,:]

    ctr = ceil(Int, (dims + 0.5)/2)

    angle = broadcast(atan2, (((1:dims[1]) - ctr[1]) ./ (dims[1]/2)), (((1:dims[2]) - ctr[2]) ./ (dims[2]/2))')
    log_rad = broadcast((x,y) -> log2(sqrt(x.^2 + y.^2)), (((1:dims[1]) - ctr[1]) ./ (dims[1]/2)), (((1:dims[2]) - ctr[2]) ./ (dims[2]/2))')
    log_rad[ctr[1], ctr[2]] = log_rad[ctr[1], ctr[2]-1]

    Xrcos, Yrcos = raisedcosine(twidth, (-twidth/2), [0 1])
    Yrcos = sqrt(Yrcos)
    YIrcos = sqrt(abs(1.0 - Yrcos.^2))

    if size(pind,1) == 2
        if any(levs .== 1)
            resdft = fftshift(fft(pyramid_subband(pyr, pind, 2)))
        else
            resdft = zeros(pind[2,:])
        end
    else
        subpyr = convert(Int, 1+prod(pind[1,:])):convert(Int, size(pyr, 1))
        resdft = reconstruct_steerable_pyramid_level(pyr[subpyr], pind[2:size(pind,1),:], log_rad, Xrcos, Yrcos, angle, nbands, levs, bands, scale)
    end

    YIrcosinterpolant = interpolate((Xrcos,), YIrcos, Gridded(Linear()))
    lo0mask = reshape(YIrcosinterpolant[log_rad], size(log_rad))
    resdft = resdft .* lo0mask

    if any(levs .== 0)
        Yrcosinterpolant = interpolate((Xrcos,), Yrcos, Gridded(Linear()))
        hi0mask = reshape(Yrcosinterpolant[log_rad], size(log_rad))
        hidft = fftshift(fft(sub_matrix(pyr, pind[1,:])))
        resdft = resdft + hidft .* hi0mask
    end

    res = real(ifft(ifftshift(resdft)))
    return res
end

# make the complex steerable pyramid real and analytic
function convert_complex_steerable_pyramid_to_real(pyr, pind; levs="all", bands="all", twidth=1, scale=0.5)
    num_scales = round(Int, log2(pind[1,1]/pind[end,1])/log2(1/scale))
    num_orientations = (size(pind,1)-2) ÷ num_scales # ÷ performs integer division in Julia

    pyrTmp = []
    ind1 = pyramid_subband_index(pind, 1)

    nband = 0
    
    for nsc in 1:num_scales
        firstBnum = (nsc-1)*num_orientations + 2

        dims = pind[firstBnum, :]
        ctr = ceil(Int, (dims+0.5)/2)
        # this is probably slow?
        ang = make_angle_grid(dims, 0, ctr)
        ang[ctr[1], ctr[2]] = -pi/2

        for nor = 1:num_orientations
            nband = (nsc-1)*num_orientations+nor+1
            ch = pyramid_subband(pyr, pind, nband)

            ang0 = pi*(nor-1)/num_orientations
            xang = mod(ang-ang0+pi, 2*pi) - pi

            # this creates an angular mask
            amask = 2*(abs(xang) .< pi/2) + (abs(xang) .== pi/2)
            amask[ctr[1], ctr[2]] = 1
            amask[:,1] = 1
            amask[1,:] = 1
            
            # and masks the fft by it
            amask = fftshift(amask)
            ch = ifft(amask.*fft(ch))
            f = 1
            ch = f * 0.5 * real(ch)

            # then creates a new pyramid
            pyrTmp = [pyrTmp; ch[:]]
        end
    end

    ind2 = pyramid_subband_index(pind, nband+1)
    pyr = [pyr[ind1]; pyrTmp; pyr[ind2]]
    
    # and returns it as a real
    return convert(Array{Float64}, real(pyr))
end

function reconstruct_complex_steerable_pyramid(pyr, pind; levs="all", bands="all", twidth=1, scale=0.5)
    pyr = convert_complex_steerable_pyramid_to_real(pyr, pind, levs=levs, bands=bands, twidth=twidth, scale=scale)
    res = reconstruct_steerable_pyramid(pyr, pind, levs=levs, bands=bands, twidth=twidth, scale=scale)
    return res
end

end
