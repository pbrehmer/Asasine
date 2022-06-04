using Revise
using Makie
using Asasine
set_theme!(ks_theme)

# KS parameters
Lx = 128
Nx = 512
dt = 1 / 16
x = Lx / Nx * collect(1:Nx)
u0 = cos.(x) .+ 0.1 * cos.(x / 16) .* (1 .+ 2 * sin.(x / 16));
ks = KSIntegrator(u0, Lx, dt);

# audio parameters
sample_rate = Float64(44100)
fps = 30
freq_step = 16
freq_idx = 10:freq_step:Nx-10
freq_func(x) = 40.0 + tanh(0.2 * x / Nx) * 4000.0
audiogen = AudioGen(sample_rate, fps, freq_func; freq_idx=freq_idx);
println("Number of pixels on y-axis: ", Nx)
println("Number of y-axis-mapped sine frequencies: ", length(freq_idx))

# start time stepping
Nt = 600
nplot = 1
att = 0.8 # additional attenuation
amp_mod = x -> x^5
start_stream(sample_rate, ks, audiogen, Nt, nplot, att, amp_mod)
