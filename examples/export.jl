using Revise
using SampledSignals
using Asasine
using Dates: now

## KS parameters
Lx = 128
Nx = 10Lx
dt = 1 / 32
x  = Lx / Nx * collect(1:Nx)
# u0 = @. 2sin(2π * x / Lx) * (cos(x) + 0.1 * cos(x / 16) * (1 + 2 * sin(x / 16)))
u0 = 10randn(Nx) .* @.(cos(x) + 0.2 * cos(x / 8) * (1 + 2 * sin(x / 4)))
# u0 = 3 * (rand(Nx) .- 0.3)
# div = 2; u0 = vcat([fill((-1.0)^n, Nx ÷ div) for n in 1:div]...)
ks = KSIntegrator{Float32}(u0, Lx, dt)

# Export static heatmap of entire KS solution
renderimage(ks, "examples/ks_L$(Lx)_$(now()).png"; Nt=1600, resolution=(1200, 800));

# Export animated heatmap
rendervideo(ks, "examples/ks_L$(Lx)_$(now()).mp4";
            Nt=2000, Ntinit=200, scroll=false, fps=60)

# Audio parameters
fps     = 12Hz
freqidx = 4:Int(Nx)÷100:Int(Nx)-4
freqs   = @. 60.0 + 3000.0tanh(0.2 * freqidx / Nx)
son     = Sonifier{Float32}(freqs, freqidx, x -> x^5)
sbuf    = SampleBuf(Float32, 44100, 1 / fps, 2);

# Export audio
sbuftot = renderaudio(sbuf, son, ks, "examples/ks_L$(Lx)_$(now()).wav");
