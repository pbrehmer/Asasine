using Revise
using SampledSignals
using Asasine

# KS parameters
Lx = 64
Nx = 4Lx
dt = 1 / 16
x  = Lx / Nx * collect(1:Nx)
u0 = @. cos(x) + 0.1 * cos(x / 16) * (1 + 2 * sin(x / 16))
ks = KSIntegrator{Float32}(u0, Lx, dt)

# Audio parameters
fps     = 30Hz
freqidx = 10:16:Int(Nx)-10
freqs   = @. 40.0 + 4000.0tanh(0.2 * freqidx / Nx)
son     = Sonifier{Float32}(freqs, freqidx, x -> x^5)
sbuf    = SampleBuf(Float32, 44100, 1 / fps, 2);

# Start live stepping
stream(sbuf, son, ks; Nt=600, att=0.8)