module Model
export make_model_data!

#using Plots ###

function calculate_data(scale::Float64, r::Float64)
	
	q = 10 .^ (-1:0.01:1) # q-values [nm^-1]
	
	#scale = 0.02 # Number of spheres per nm^3 [nm^-3]
	#r = 0.5 # nm
	rho = 1 # nm^-2
	background = 0 # nm^-1
	
	qr = q .* r
	intensity = scale .* ((3*(4*pi*r^3/3)*rho).*(sin.(qr) - qr .* cos.(qr)) ./ qr.^3).^2 .+ background # nm^-1
	
	phi = scale * (4*pi*r^3/3)
	
	psi = 3*phi / (1 - phi)
	b = psi .* ((cos.(qr) + qr .* sin.(qr)) .* (sin.(qr) - qr .* cos.(qr))) ./ qr.^3
	c = psi .* (sin.(qr) - qr .* cos.(qr)).^2 ./ qr.^3
	d = 1 .+ psi .* (qr.^2 .* sin.(qr) .* cos.(qr)) ./ qr.^3
	e = psi .* (qr.^2 .* (sin.(qr).^2)) ./ qr.^3
	f = psi .* (qr .* sin.(qr) .* (sin.(qr) - qr .* cos.(qr))) ./ qr.^3
	g = - psi .* (qr .* cos.(qr) .* (sin.(qr) - qr .* cos.(qr))) ./ qr.^3
	
	X = 1 .+ b + (2 .* e .* f .* g + d .* (f.^2 - g.^2)) ./ (d.^2 + e.^2)
	Y = c + (2 .* d .* f .* g - e .* (f.^2 - g.^2)) ./ (d.^2 + e.^2)
	
	S = (Y ./ c) ./ (X.^2 + Y.^2)
	
	#gr() ###
	#global plt = plot(q, intensity; xaxis = :log, yaxis = :log) ###
	
	intensity = S .* intensity
	
	#plot!(q, intensity)
	
	return q::Array{Float64, 1}, intensity::Array{Float64, 1}
	
end

function save_data_to_file!(data::Tuple{Array{Float64, 1}, Array{Float64, 1}}, outputFileName::String)
	
	io = open(outputFileName, "w")
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

function make_model_data!(scale::Float64, r::Float64)
	
	save_data_to_file!(calculate_data(scale::Float64, r::Float64), "calculated.data")
	
	return nothing
	
end

end # Model
