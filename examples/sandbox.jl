using Asasine

# KS parameters
Lx = 80
Nx = 1024
dt = 1 / 16
x = Lx / Nx * collect(1:Nx)
u0 = cos.(x) .+ 0.1 * cos.(x / 16) .* (1 .+ 2 * sin.(x / 16));
# u0 = 3.0 * exp.(-(x .- Nx/2).^2 ./ sqrt(Nx))
# u0 = 2.0 * (rand(length(x)) .- 0.5)
# u0 = 2.0 * rand(length(x))
ks = KSIntegrator(u0, Lx, dt);

# audio parameters
sample_rate = Float64(44100)
fps = 30
downsampler = 16
println("Video resolution per time step: ", Nx)
println("Audio resolution per time step: ", Nx / downsampler)
# freqs = collect(range(40, 400, Nx))
freqs = [50.0 + tanh(0.3 * i / Nx) * 2500.0 for i = 0:downsampler:Nx-1]
audiogen = AudioGen(sample_rate, fps, freqs);

# start time stepping
Nt = 800
nplot = 1
att = 0.6 # additional attenuation
amp_mod = x -> x^5
start(sample_rate, ks, audiogen, Nt, nplot, att, amp_mod)
