module SampleMaker

##########################################################################

module Tools
export make_hard_sphere_samples, set_sld_offsets, sld_scan
using FileIO, ColorTypes, Serialization

global const PX_PER_NM = 10
global const SAMPLE_RADIUS_NM = 10

global const SAMPLE_RADIUS_PX = PX_PER_NM * SAMPLE_RADIUS_NM

function poisson(lambda::Float64)
	
	k = 1
	p = rand()
	q = exp(-lambda)
	
	while p >= q
		p *= rand()
		k += 1
	end
	
	return k::Int
	
end

function empty_sample()
	
	sample = -Inf * ones(2*SAMPLE_RADIUS_PX+1, 2*SAMPLE_RADIUS_PX+1, 2*SAMPLE_RADIUS_PX+1)
	center = (SAMPLE_RADIUS_PX+1, SAMPLE_RADIUS_PX+1, SAMPLE_RADIUS_PX+1)
	
	x1 = repeat(1:2*SAMPLE_RADIUS_PX+1, 1, 2*SAMPLE_RADIUS_PX+1, 2*SAMPLE_RADIUS_PX+1)
	x2 = permutedims(x1, [2 1 3])
	x3 = permutedims(x1, [3 2 1])
	
	dists = sqrt.( (x1.-center[1]).^2 + (x2.-center[2]).^2 + (x3.-center[3]).^2 )
	sample[abs.(dists) .<= SAMPLE_RADIUS_PX] .= 0
	
	sample = Float16.(sample)
	
	return sample::Array{Float16, 3}
	
end

function sphere(center::Tuple{Int, Int, Int}, radius::Float64)
	
	sphere = zeros(2*SAMPLE_RADIUS_PX+1, 2*SAMPLE_RADIUS_PX+1, 2*SAMPLE_RADIUS_PX+1)
	
	x1 = repeat(1:2*SAMPLE_RADIUS_PX+1, 1, 2*SAMPLE_RADIUS_PX+1, 2*SAMPLE_RADIUS_PX+1)
	x2 = permutedims(x1, [2 1 3])
	x3 = permutedims(x1, [3 2 1])
	
	dists = sqrt.( (x1.-center[1]).^2 + (x2.-center[2]).^2 + (x3.-center[3]).^2 )
	sphere[abs.(dists) .<= radius] .= 1
	
	sphere = Float16.(sphere)
	
	return sphere::Array{Float16, 3}
	
end

function sampRms_for_offset(samplePath::String, offset::Float64)
	
	sampRms = 0
	
	n = 1
	while isfile(samplePath*string(n, pad = 5)*".samp")
		
		sample = deserialize(samplePath*string(n, pad = 5)*".samp")
		sample = sample .+ Float16(offset)
		
		let temp = Float64.(deepcopy(sample))
			temp[temp .== -Inf] .= 0
			sampRms += (sum(temp) * (1/PX_PER_NM)^3)^2
		end
		
		n += 1
		
	end
	n -= 1
	
	sampRms /= n
	sampRms = sqrt(sampRms)
	
	return sampRms::Float64
	
end

function set_scale()
	io = open("sample/scale.meta", "w")
	println(io, "Pixels per nanometer: ", PX_PER_NM)
	println(io, "Sample radius in nanometers: ", SAMPLE_RADIUS_NM)
	close(io)
	return nothing
end

function set_sld_offsets(targetIntensity::Float64)
	
	println("\nSetting SLD offsets...")
	
	io = open("sample/offsets.meta", "w")
	println(io, "SLD offset [nm^-2]")
	close(io)
	
	samplePath = "sample/"
	
	targetIntensity *= (4 * pi * SAMPLE_RADIUS_NM^3 / 3)
	
	n = 1
	while isfile(samplePath*string(n, pad = 5)*".samp")
		
		sample = deserialize(samplePath*string(n, pad = 5)*".samp")
		
		sampleSldSum = 0
		let temp = Float64.(deepcopy(sample))
			temp[temp .== -Inf] .= 0
			sampleSldSum = sum(temp) * (1/PX_PER_NM)^3
		end
		
		offset = - (sampleSldSum - sqrt(targetIntensity)) / (4 * pi * SAMPLE_RADIUS_NM^3 / 3)
		
		io = open("sample/offsets.meta", "a")
		println(io, offset)
		close(io)
		
		n += 1
		
	end
	
	println("...done.")
	
	return nothing
	
end

function offset_scan(offsets::AbstractRange{Float64})
	println("\nSLD offset [nm^-2], Intensity at q = 0 [nm^2]")
	intensities = zeros(0)
	for offset in offsets
		intensity = sampRms_for_offset("sample/", offset)^2
		intensities = [intensities; [intensity]]
		println(offset, ", ", intensity)
	end
	return intensities::Array{Float64, 1}
end

function write_sample_files(sample::Array{Float16, 3}, id::String)
	
	function extract_image_from_file(fileName::String)
		if !isfile(fileName)
			println("""The file \""""*fileName*"""\" is missing.""")
			exit()
		else
			image = RGBA.(load(fileName))
			(r, g, b, a) = (red.(image), green.(image), blue.(image), alpha.(image))
			image = cat(r, g, b, a; dims = 3)
			image = round.(UInt8, 255*image)
			return image::Array{UInt8, 3}
		end
	end
	
	function make_figure(sample::Array{Float16, 3})
		
		function project(sample::Array{Float16, 3}, dimension::Int)
			
			projection = Float64.(deepcopy(sample))
			projection[projection .< 0] .= 0
			projection = sum(projection; dims = dimension)
			projection = dropdims(projection; dims = dimension)
			projection[projection .< 0] .= 0
			
			return projection::Array{Float64, 2}
			
		end
		
		mask = sample[:, :, Int((size(sample)[3] - 1)/2 + 1)]
		mask[mask .>= 0] .= 0.0
		mask[mask .< 0] .= 1.0
		mask = [mask mask mask]
		
		x1Projection = project(sample, 1)
		x2Projection = project(sample, 2)
		x3Projection = project(sample, 3)
		
		annotation = zeros(size(mask))
		annotation[1:25, 1:25] = Float64.(extract_image_from_file("../rsc/1.png")[:, :, 1]) ./ 255
		annotation[1:25, 202:226] = Float64.(extract_image_from_file("../rsc/2.png")[:, :, 1]) ./ 255
		annotation[1:25, 403:427] = Float64.(extract_image_from_file("../rsc/3.png")[:, :, 1]) ./ 255
		annotation[191:196, 579:598] .= 1
		
		figure = [x3Projection x2Projection x1Projection]
		if maximum(figure) != 0
			figure = figure ./ maximum(figure)
		end
		figure = 1 .- figure
		
		figure = RGB.(figure - 0.1 .* mask - 0.8 .* annotation, figure, figure)
		
		return figure::Array{RGB{Float64}, 2}
		
	end
	
	figure = make_figure(sample)
	save("sample/"*id*".png", figure)
	save("../../output/sample_projections/"*id*".png", figure)
	io = open("sample/"*id*".samp", "w")
	serialize(io, sample)
	close(io)
	
end

function make_hard_sphere_samples(numOfSamples::Int, numOfSpheresPerNm3::Float64, sphereRadiusNm::Float64)
	
	rm("sample/"; recursive = true, force = true)
	mkdir("sample/")
	
	sphereRadius = PX_PER_NM * sphereRadiusNm
	
	boxSideLength = 2*(SAMPLE_RADIUS_PX + sphereRadius + 1) # px
	boxVolume = boxSideLength^3 # px^3
	
	boxVolumeNm3 = boxVolume * (1/PX_PER_NM)^3 # nm^3
	
	numOfSpheres = boxVolumeNm3 * numOfSpheresPerNm3
	
	set_scale()
	
	io = open("sample/description.txt", "w")
	print(io, "The sample is made of uniform hard spheres with a radius of ", round(sphereRadiusNm; sigdigits = 5), " nm.")
	print(io, " The concentration of spheres is ", round(numOfSpheresPerNm3; sigdigits = 5), " per nm^3.")
	print(io, " This corresponds to an average of ", round(numOfSpheresPerNm3*(4 * pi * SAMPLE_RADIUS_NM^3 / 3); sigdigits = 5), " spheres per sample volume,")
	print(io, " a volume fraction of ", round(numOfSpheresPerNm3*(4 * pi * sphereRadiusNm^3 / 3); sigdigits = 5), ",")
	print(io, " and a molar concentration of ", round((10^3)*(10.0^24)*numOfSpheresPerNm3/(6.02214076*10.0^23); sigdigits = 5), " mM.")
	print(io, " The scale bar shows ", round(20/PX_PER_NM; sigdigits = 5), " nm.")
	close(io)
	mkdir("../../output/sample_projections/")
	cp("sample/description.txt", "../../output/sample_projections/description.txt")
	
	println()
	
	for n in 1:numOfSamples
		
		println("Making sample ", n, " of ", numOfSamples, "...")
		
		sample = empty_sample()
		
		midpoint = 101
		sideRange = midpoint-Int(boxSideLength/2):midpoint+Int(boxSideLength/2)
		for i in 1:poisson(numOfSpheres)
			sphereCenter = (rand(sideRange), rand(sideRange), rand(sideRange))
			newSample = sample + sphere(sphereCenter, sphereRadius)
			while maximum(newSample) > 1
				sphereCenter = (rand(sideRange), rand(sideRange), rand(sideRange))
				dist = sqrt.( (midpoint-sphereCenter[1]).^2 + (midpoint-sphereCenter[2]).^2 + (midpoint-sphereCenter[3]).^2 )
				while dist > SAMPLE_RADIUS_PX - sphereRadius
					sphereCenter = (rand(sideRange), rand(sideRange), rand(sideRange))
					dist = sqrt.( (midpoint-sphereCenter[1]).^2 + (midpoint-sphereCenter[2]).^2 + (midpoint-sphereCenter[3]).^2 )
				end
				newSample = sample + sphere(sphereCenter, sphereRadius)
			end
			sample = deepcopy(newSample)
		end
		
		write_sample_files(sample, string(n, pad = 5))
		
	end
	
	println("Done.")
	
	return nothing
	
end

end # Tools

##########################################################################

end # SampleMaker
