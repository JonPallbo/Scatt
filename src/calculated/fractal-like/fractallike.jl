function riemann_zeta(s::Float64)
	
	degree = 1E6
	z = 0
	for n in 1:degree
		z += n^-s
	end
	
	return z::Float64
	
end

function object(d_object::Float64, q::Array{Float64, 1}, r::Float64)
	
	c = r^(2*(3-d_object))
	
	intensity = min.(c .* (r.*q).^-(d_object), c .* (r.*q).^-4)
	
	return intensity::Array{Float64, 1}
	
end

function fractal(d_fractal::Float64, d_object::Float64, r_monomer::Float64, form_factor)
	
	intensity = 0 .* form_factor(1.0)
	p = (d_fractal+3*d_object-9)/(3-d_object)
	z = riemann_zeta(-(p+1-(4/(3-d_object)-1)))
	for n in 1:1E6
		r = r_monomer * n^(1/(3-d_object))
		N = (n^p)/z
		intensity .+= N .* form_factor(r)
	end
	return intensity::Array{Float64, 1}
	
end

##########################################################################

using Plots

q = 10 .^ (-1:0.1:1) # [nm^-1]

gr()
plt = plot()
plot!(q[1:6], 1E2.*q[1:6].^-2.5;
	seriestype = :line,
	xaxis = :log,
	yaxis = :log,
	minorgrid = true,
	xticks = 10.0 .^ (-1:1),
	yticks = 10.0 .^ (-5:6),
	linewidth = 1.0,
	xlabel = "Q-value [nm\u207B\u00B9]",
	ylabel = "Intensity [a.u.]",
	foreground_color_legend = nothing,
	background_color_legend = nothing,
	framestyle = :box,
	legend = :topright,
	legendfontsize = 7,
	guidefontsize = 7,
	tickfontsize = 7,
	linestyle = :solid,
	linecolor = RGB(0, 1, 0),
	linealpha = 1.0,
	lab = "\u221D q\u207B\u00B2\u22C5\u2075")

plot!([NaN], [NaN];
	linealpha = 0,
	lab = " ")

plot!([NaN], [NaN];
	linestyle = :dot,
	linecolor = RGB(0, 0, 0),
	linealpha = 0.5,
	lab = "Normal")

plot!([NaN], [NaN];
	linestyle = :solid,
	linecolor = RGB(0, 0, 0),
	linealpha = 0.5,
	lab = "Fractal-like")

plot!([NaN], [NaN];
	linealpha = 0,
	lab = "(d = 2.5)")

plot!([NaN], [NaN];
	linealpha = 0,
	lab = " ")

plot!([1], [1];
	seriestype = :scatter,
	markershape = :circle,
	markersize = 1,
	markercolor = RGB(1, 0, 1),
	markerstrokecolor = RGB(1, 0, 1),
	markeralpha = 0.5,
	lab = "Cylindrish")

plot!([1], [1];
	seriestype = :scatter,
	markershape = :circle,
	markersize = 1,
	markercolor = RGB(0, 0, 1),
	markerstrokecolor = RGB(0, 0, 1),
	markeralpha = 0.5,
	lab = "Spherish")

plot!(q, 1E1.*object(1.0, q, 1.0);
	linestyle = :dot,
	linecolor = RGB(1, 0, 1),
	linealpha = 0.5,
	lab = nothing)

plot!(q, 1E0.*object(0.0, q, 1.0);
	linestyle = :dot,
	linecolor = RGB(0, 0, 1),
	linealpha = 0.5,
	lab = nothing)

plot!(q, 1E1.*fractal(2.5, 1.0, 1.0, r -> object(1.0, q, r));
	linestyle = :solid,
	linecolor = RGB(1, 0, 1),
	linealpha = 0.5,
	lab = nothing)

plot!(q, 1E0.*fractal(2.5, 0.0, 1.0, r -> object(0.0, q, r));
	linestyle = :solid,
	linecolor = RGB(0, 0, 1),
	linealpha = 0.5,
	lab = nothing)

savefig(plt, "figure.pdf")
