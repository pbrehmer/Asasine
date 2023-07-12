using Revise
using SampledSignals
using Asasine
using Dates: now

# KS parameters
Lx = 128
Nx = 8Lx
dt = 1 / 32
x  = Lx / Nx * collect(1:Nx)
# u0 = @. 2sin(2π * x / Lx) * (cos(x) + 0.1 * cos(x / 16) * (1 + 2 * sin(x / 16)))
u0 = 5randn(Nx) .* @.(cos(x) + 0.1 * cos(x / 16) * (1 + 2 * sin(x / 16)))
# u0 = 4rand(Nx)
ks = KSIntegrator{Float32}(u0, Lx, dt)

# Audio parameters
fps     = 30Hz
freqidx = 10:16:Int(Nx)-10
freqs   = @. 40.0 + 4000.0tanh(0.2 * freqidx / Nx)
son     = Sonifier{Float32}(freqs, freqidx, x -> x^5)
sbuf    = SampleBuf(Float32, 44100, 1 / fps, 2)

# TODO: Export audio

# Export static heatmap of entire KS solution
renderimage(ks, "examples/ks_L$(Lx)_$(now()).png"; Nt=1400, resolution=(1200, 800))

# TODO: Export animated heatmap
