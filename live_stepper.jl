import FFTW, GLMakie, PortAudio
using Makie: PriorityObservable
GLMakie.set_theme!(GLMakie.theme_black())

Base.@kwdef mutable struct KSIntegrator
    u0::Vector{Float64}
    Lx::Float64
    dt::Float64

    ### additional physical parameters
    Nx::Integer = length(u0)
    kx::Vector{Integer} = vcat(0:Nx/2-1, 0, -Nx/2+1:-1) # shifted PBC integer wavenumbers
    alpha::Vector{Float64} = 2pi/Lx * kx

    # convenience variables
    L::Vector{Float64} = alpha.^2 - alpha.^4
    G::Vector{ComplexF64} = -0.5im * alpha
    dt2::Float64 = dt/2
    dt32::Float64 = 3*dt/2
    A::Vector{Float64} = ones(Nx) + dt2*L
    B::Vector{Float64} = (ones(Nx) - dt2*L).^-1
    
    # compute in-place planned FFTs
    FFT! = FFTW.plan_fft!((1+0im)*u0, flags=FFTW.ESTIMATE)
    IFFT! = FFTW.plan_ifft!((1+0im)*u0, flags=FFTW.ESTIMATE)

    # prepare time evolution in Fourier domain
    Nn::Vector{ComplexF64} = G .* FFTW.fft(u0.^2)
    Nn1::Vector{ComplexF64} = copy(Nn)
    u::Vector{ComplexF64} = FFT! * u0
end

Base.@kwdef mutable struct AudioGen
        sample_rate::Float64
        fps::Integer
        freqs::Vector{Float64}
        last_phases = zeros(Float64, length(freqs))
        grain_time::Float64 = Float64(1 / fps)
        buf_size::Integer = ceil(Integer, grain_time * sample_rate)
        buf::Matrix{Float64} = zeros(Float64, buf_size, 2)
end

include("./src/KS.jl")
include("./src/Audio.jl")

### main ###
# KS parameters
Lx = 64
Nx = 100
dt = 1/16
x = Lx/Nx * collect(1:Nx)
u0 = cos.(x) .+ 0.1*cos.(x/16) .* (1 .+ 2*sin.(x/16));
# u0 = 3.0 * exp.(-(x .- Nx/2).^2 ./ sqrt(Nx))
# u0 = 2.0 * (rand(length(x)) .- 0.5)
# u0 = 2.0 * rand(length(x))
ks = KSIntegrator(u0 = u0, Lx = Lx, dt = dt)

# audio parameters 
sample_rate = Float64(44100)
fps = 30
# freqs = collect(range(40, 400, Nx))
freqs = [70.0 + tanh(0.35*i/Nx) * 7000.0 for i in 0:Nx-1]
audiogen = AudioGen(sample_rate = sample_rate, fps = fps, freqs = freqs)

# start time stepping
Nt = 800
nplot = 1
att = 0.6 # additional attenuation
amp_mod = x -> x^5
start(ks, audiogen, Nt, nplot, att, amp_mod)
