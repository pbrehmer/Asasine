using Revise
using Makie
using Asasine
set_theme!(theme_black())

# KS parameters
Lx = 80
Nx = 512
dt = 1 / 32
x = Lx / Nx * collect(1:Nx)
u0 = cos.(x) .+ 0.1 * cos.(x / 16) .* (1 .+ 2 * sin.(x / 16));
# u0 =
#     4.0 * exp.(-(x .- 3 * x[Nx÷2] / 4) .^ 2 / sqrt(Nx)) -
#     2.0 * exp.(-(x .- x[Nx÷2] / 4) .^ 2 / sqrt(Nx))
# u0 = 2.0 * (rand(length(x)) .- 0.5)
ks = KSIntegrator(u0, Lx, dt);

# audio parameters
sample_rate = Float64(44100)
fps = 30
freq_step = 16
freq_idx = 10:freq_step:Nx-10
freqs = [50.0 + tanh(0.3 * i / Nx) * 10000.0 for i in freq_idx]
audiogen = AudioGen(sample_rate, fps, freqs; freq_idx=freq_idx);
println("Number of pixels on y-axis: ", Nx)
println("Number of y-axis-mapped sine frequencies: ", length(freq_idx))

# start time stepping
Nt = 600
nplot = 1
att = 0.6 # additional attenuation
amp_mod = x -> x^5
start_stream(sample_rate, ks, audiogen, Nt, nplot, att, amp_mod)
