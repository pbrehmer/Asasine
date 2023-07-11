using Revise
using GLMakie
using SampledSignals
using Asasine
set_theme!(kstheme)
GLMakie.activate!(inline=false)  # Open visuals in window

# KS parameters
Lx = Float32(64)
Nx = 4Lx
dt = Float32(1 / 16)
x  = Lx / Nx * collect(1:Nx)
u0 = Float32.(@. cos(x) + 0.1 * cos(x / 16) * (1 + 2 * sin(x / 16)))

# Audio parameters
fps     = 30Hz
freqidx = 10:16:Int(Nx)-10
freqs   = Float32.(@. 40.0 + 4000.0tanh(0.2 * freqidx / Nx))
son     = Sonifier(freqs, freqidx, x -> x^5)
println("Number of pixels on y-axis: ", Nx)
println("Number of y-axis-mapped sine frequencies: ", length(freqidx))

## Initialize integrator and buffer, and start live stepping
ks = KSIntegrator{Float32,ComplexF32}(; u0, Lx, dt);
sbuf = SampleBuf(Float32, 44100, 1 / fps, 2)
stream(sbuf, son, ks; Nt=600, att=0.8)