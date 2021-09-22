rm("../output/"; recursive = true, force = true)
mkdir("../output/")
include("calculated/model.jl")
include("simulated/sampleMaker.jl")
include("scatt.jl")

cd("calculated/")

using ..Model
make_model_data!(0.005, 1.5)

cd("..")

using ..Scatt.Finitizer
finitize_data("calculated/calculated.data")
targetIntensity = intensity_at_zero("../output/calculated_finitized.data")

cd("simulated/")

using ..SampleMaker.Tools
make_hard_sphere_samples(100, 0.005, 1.5)
set_sld_offsets!(targetIntensity)

cd("..")

using ..Scatt.Scatterer
measure!("simulated/sample/")

using ..Scatt.Plotter
make_plot!()
