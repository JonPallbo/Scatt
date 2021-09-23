module Scatt

##########################################################################

module Auxiliaries
export read_data_from_file, save_data_to_file

function read_data_from_file(FileName::String)
	
	dataString = read(FileName, String)
	
	dataString = replace(dataString, "\r\n" => "\n")
	dataString = replace(dataString, "\r" => "\n")
	dataString = replace(dataString, "\t" => "")
	dataString = replace(dataString, " " => "")
	dataString = String.(split(dataString, "\n"; keepempty = false))
	dataString = dataString[2:end]
	
	qValues = zeros(length(dataString))
	intensities = zeros(length(qValues))
	
	for i in 1:length(dataString)
		s = dataString[i]
		qValues[i] = parse(Float64, split(s, ",")[1])
		intensities[i] = parse(Float64, split(s, ",")[2])
	end
	
	return qValues::Array{Float64, 1}, intensities::Array{Float64, 1}
	
end

function save_data_to_file(data::Tuple{Array{Float64, 1}, Array{Float64, 1}}, outputFileName::String)
	
	io = open(outputFileName, "a")
	if isempty(read(outputFileName, String))
		println(io, "q-Value [nm^-1], Intensity [nm^-1]")
	end
	
	for i in 1:length(data[1])
		
		local qValue = data[1][i]
		local intensity = data[2][i]
		
		println(io, qValue, ", ", intensity)
		
	end
	
	close(io)
	
	return nothing
	
end

end # Auxiliaries

##########################################################################

module Finitizer
export finitize_data, intensity_at_zero
using ..Auxiliaries, FFTW, Interpolations

function expand_data(data::Tuple{Array{Float64, 1}, Array{Float64, 1}}, qValues::Array{Float64, 3})
	
	intensity = zeros(size(qValues))
	intensity[:] = LinearInterpolation(data[1], data[2]; extrapolation_bc = ((Line(), Flat()),))[qValues[:]]
	
	return intensity::Array{Float64, 3}
	
end

function finitize_data(filePath::String, qRange::Array{Float64, 1}, sampleRadiusNm::Float64)
	
	data = read_data_from_file(filePath)
	(qs, ints) = finitize_data(data, qRange, sampleRadiusNm)
	rm("../output/finitized.data"; force = true)
	save_data_to_file((qs, ints), "../output/finitized.data")
	
	return nothing
	
end

function finitize_data(data::Tuple{Array{Float64, 1}, Array{Float64, 1}}, qs::Array{Float64, 1}, sampleRadiusNm::Float64)
	
	println("\nFinitizing data...")
	
	function reduce_data(q::Array{Float64, 3}, data::Array{Float64, 3}, px_per_inv_nm::Int, qRange::Array{Float64, 1})
		
		sideLength = size(data, 1)
		midpoint = Int((sideLength-1)/2 + 1)
		
		lowQ = [q[midpoint:midpoint+10, midpoint:midpoint+10, midpoint:midpoint+10][:] data[midpoint:midpoint+10, midpoint:midpoint+10, midpoint:midpoint+10][:]]
		lowQ = sortslices(lowQ; dims = 1)
		lowQ = lowQ[indexin(unique(lowQ[:, 1]), lowQ[:, 1]), :]
		lowQ = lowQ[lowQ[:, 1] .< 11/px_per_inv_nm, :]
		rs = [lowQ[:, 1]; ((11:Int((sideLength-1)/2)) ./ px_per_inv_nm)]
		intensities = [lowQ[:, 2]; data[midpoint+11:end, midpoint, midpoint]]
		intensities = LinearInterpolation(rs, intensities)[qRange]
		
		return intensities::Array{Float64, 1}
		
	end
	
	scale = 1
	r = sampleRadiusNm
	background = 0
	
	px_per_inv_nm = round(Int, 2*r)
	
	sideLength = 401
	x1 = repeat(1:sideLength, 1, sideLength, sideLength)
	x2 = permutedims(x1, [2 1 3])
	x3 = permutedims(x1, [3 2 1])
	midpoint = Int((sideLength-1)/2 + 1)
	q = sqrt.( (x1.-midpoint).^2 + (x2.-midpoint).^2 + (x3.-midpoint).^2 ) ./ px_per_inv_nm # nm^-1
	
	qr = q .* r
	mask = scale .* ((3*1).*(sin.(qr) - qr .* cos.(qr)) ./ qr.^3).^2 .+ background # nm^2
	mask[midpoint, midpoint, midpoint] = 1
	
	mask[abs.(q .* r) .> 98.950063] .= 0
	
	mask = mask ./ sum(mask)
	
	intensity = expand_data(data, q)
	
	convolution = abs.(ifftshift(ifft(fft(mask) .* fft(intensity))))
	
	ints = reduce_data(q, convolution, px_per_inv_nm, qs)
	
	println("...done.")
	
	return qs::Array{Float64, 1}, ints::Array{Float64, 1}
	
end

function intensity_at_zero(filePath::String)
	
	data = read_data_from_file(filePath)
	
	intensity = 0
	for i in 2:min(10, length(data[1]))
		intensity += data[2][1] + (data[2][i] - data[2][1])/(data[1][i] - data[1][1]) * (0 - data[1][1])
	end
	intensity /= min(10, length(data[1])) - 1
	
	return intensity::Float64
	
end

end # Finitizer

##########################################################################

module Scatterer
export get_sampsum_for_offset, measure
using Serialization

function load_meta(samplePath::String)
	
	let dataString = read(samplePath*"scale.meta", String)
		
		dataString = replace(dataString, "\r\n" => "\n")
		dataString = replace(dataString, "\r" => "\n")
		dataString = replace(dataString, "\t" => "")
		dataString = replace(dataString, " " => "")
		dataString = String.(split(dataString, "\n"; keepempty = false))
		
		global PX_PER_NM = parse(Float64, split(dataString[1], ":")[2])
		global SAMPLE_RADIUS_NM = parse(Float64, split(dataString[2], ":")[2])
		
		global SAMPLE_RADIUS_PX = PX_PER_NM * SAMPLE_RADIUS_NM
		
	end
	
	let dataString = read(samplePath*"offsets.meta", String)
		
		dataString = replace(dataString, "\r\n" => "\n")
		dataString = replace(dataString, "\r" => "\n")
		dataString = replace(dataString, "\t" => "")
		dataString = replace(dataString, " " => "")
		dataString = String.(split(dataString, "\n"; keepempty = false))
		
		global SLD_OFFSETS = parse.(Float64, dataString[2:end])
		
	end
	
	return nothing
	
end

function make_probes(qs::Array{Tuple{Float64, Float64, Float64}, 1})
	
	probes = Array{Array{Complex{Float16}, 3}, 1}()
	
	x1 = Float16.(repeat(1:Int(2*SAMPLE_RADIUS_PX+1), 1, Int(2*SAMPLE_RADIUS_PX+1), Int(2*SAMPLE_RADIUS_PX+1)))
	x2 = permutedims(x1, [2 1 3])
	x3 = permutedims(x1, [3 2 1])
	
	for i in 1:length(qs)
		q = qs[i]
		probe = Complex{Float16}.(exp.((x1.*q[1]+x2.*q[2]+x3.*q[3]).*im))
		probes = [probes; [probe]]
		println(i, "/", length(qs))
	end
	
	return probes::Array{Array{Complex{Float16}, 3}, 1}
	
end

function get_intensity(sample::Array{Float16, 3}, probe::Array{Complex{Float16}, 3})
	
	sample = deepcopy(sample)
	sample[sample .== -Inf] .= 0
	
	intensity = abs(sum(Float64.(sample) .* Complex{Float64}.(probe) .* (1/PX_PER_NM)^3))^2
	
	return intensity::Float64
	
end

function measure(samplePath::String)
	
	load_meta(samplePath)
	
	sampleVolume = (4 * pi * SAMPLE_RADIUS_NM^3 / 3)
	
	numMeasur = 3
	for k in 1:numMeasur
		
		println("\nSetting up measurement ", k, " of ", numMeasur, "...\n")
		
		global qs = Array{Tuple{Float64, Float64, Float64}, 1}()
		rs = 10 .^ (-1:0.01:1) # q-range [nm^-1]
		for r in rs
			local q = (1, 1, 1)
			while sqrt(sum(q.^2)) > 1
				q = (2*rand()-1, 2*rand()-1, 2*rand()-1)
			end
			q = r .* (q ./ sqrt(sum(q.^2)))
			q = q ./ PX_PER_NM # pixel wave vector
			global qs = [qs; [q]]
		end
		
		probes = make_probes(qs)
		
		println("\n...done.")
		
		intensities = zeros(length(rs))
		
		n = 1
		while isfile(samplePath*string(n, pad = 5)*".samp")
			
			sample = deserialize(samplePath*string(n, pad = 5)*".samp")
			sample = sample .+ Float16(SLD_OFFSETS[n])
			
			println("\nq-value [nm^-1], Intensity [nm^-1]", " # Measurement ", k, " of ", numMeasur, ".")
			for i in 1:length(rs)
				intensities[i] += get_intensity(sample, probes[i])
				local intensity = intensities[i]
				local qValue = sqrt(sum((qs[i].*PX_PER_NM).^2))
				println(qValue, ", ", intensity/sampleVolume, " # "*string(n, pad = 5)*" [", i, "/", length(rs), "]")
			end
			
			n += 1
			
		end
		n -= 1
		
		intensities = intensities ./ (n * sampleVolume)
		
		outputFileName = "../output/simulated.data"
		io = open(outputFileName, "a")
		if isempty(read(outputFileName, String))
			println(io, "q-Value [nm^-1], Intensity [nm^-1]")
		end
		
		for i in 1:length(rs)
			
			local intensity = intensities[i]
			local qValue = sqrt(sum((qs[i].*PX_PER_NM).^2))
			
			println(io, qValue, ", ", intensity)
			
		end
		
		close(io)
		
	end
	
	return nothing
	
end

end # Scatterer

##########################################################################

module Plotter
export make_plot
using ..Auxiliaries, ..Finitizer, Plots

function make_plot()
	
	gr()
	
	plt = plot()
	
	(qValues, intensities) = read_data_from_file("../output/simulated.data")
	plot!(qValues, intensities;
		seriestype = :scatter,
		xaxis = :log,
		yaxis = :log,
		minorgrid = true,
		markershape = :circle,
		markersize = 2,
		markercolor = RGB(0, 0, 1),
		markerstrokecolor = RGB(0, 0, 1),
		markeralpha = 0.5,
		xlabel = "q-Value [nm\u207B\u00B9]",
		ylabel = "Intensity [nm\u207B\u00B9]",
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		lab = "Simulated",
		framestyle = :box)
	
	(qValues, intensities) = read_data_from_file("../output/finitized.data")
	plot!(qValues, intensities;
		seriestype = :line,
		linewidth = 3,
		linecolor = RGB(0, 1, 0),
		linealpha = 0.5,
		lab = "Calculated")
	
	savefig(plt, "../output/figure.pdf")
	
	return nothing
	
end

end # Plotter

##########################################################################

end # Scatt
