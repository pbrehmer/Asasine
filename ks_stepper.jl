import FFTW, GLMakie, PortAudio
using Makie: PriorityObservable
GLMakie.set_theme!(GLMakie.theme_black())


Base.@kwdef mutable struct KSIntegrator
    u0::Vector{Float64}
    Lx::Float64
    dt::Float64

    ### additional physical parameters
    Nx::Int64 = length(u0)
    kx::Vector{Int64} = vcat(0:Nx/2-1, 0, -Nx/2+1:-1) # shifted PBC integer wavenumbers
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

function evolve!(ks::KSIntegrator, steps::Int=1)
    for n = 0:steps-1
        # compute Nn = G .* fft(real(ifft(u)).^2) inbounds with fastmath
        for i = 1:length(ks.Nn)
            @inbounds ks.Nn1[i] = ks.Nn[i]
            @inbounds ks.Nn[i] = ks.u[i]
        end

        ks.IFFT! * ks.Nn

        for i = 1:length(ks.Nn)
            @fastmath @inbounds ks.Nn[i] = ks.Nn[i] * ks.Nn[i]
        end
        
        ks.FFT! * ks.Nn
        
        for i = 1:length(ks.Nn)
            @fastmath @inbounds ks.Nn[i] = ks.G[i] * ks.Nn[i]
        end

        # compute u = B .* (A .* u + dt32*Nn - dt2*Nn1) inbounds with fastmath
        for i = 1:length(ks.u)
            @fastmath @inbounds ks.u[i] = ks.B[i] *
                (ks.A[i] * ks.u[i] + ks.dt32*ks.Nn[i] - ks.dt2*ks.Nn1[i])
        end
    end
end

function getSolution(ks::KSIntegrator)
    return real(FFTW.ifft(ks.u))
end

function ksIntegrate!(ks::KSIntegrator, Nt::Int64, nplot::Int64)
    nsave = round(Int64, Nt/nplot)+1 # total number of saved time steps
    x = ks.Lx/ks.Nx * collect(1:ks.Nx)
    t = ks.dt*nplot * collect(0:nsave-1)
    U = zeros(Float64, nsave, ks.Nx)
    U[1,:] = ks.u0

    evolve!(ks)
    for i in 2:nsave
        U[i, :] = getSolution(ks)
        evolve!(ks, nplot)
    end

    U, t, x
end

function getSineStack(u::Vector{Float64}, freqs::Vector{Float64}, amp_mapping::Function,
        last_phases::Vector{Float64}, sample_rate::Real, grain_time::Real)
    # mono version
    buf_size = ceil(Int64, grain_time * sample_rate)
    buf = zeros(Float64, buf_size)
    for i = 1:length(freqs)
        buf = buf .+ amp_mapping(u[i]) .*
            sin.(last_phases[i] .+ 2pi * (0:buf_size-1) * freqs[i] / sample_rate)
        last_phases[i] = (last_phases[i] + 2pi * buf_size * freqs[i] / sample_rate) % 2pi
    end
    buf, last_phases
end

function getSineStack(u::Matrix{Float64}, freqs::Vector{Float64}, amp_mapping::Function,
        last_phases::Vector{Float64}, sample_rate::Real, grain_time::Real)
    # stereo version
    buf_size = ceil(Int64, grain_time * sample_rate)
    buf = zeros(Float64, buf_size, 2)

    for i = 1:length(freqs)
        buf[:, 1] = buf[:, 1] .+ amp_mapping(u[i, 1]) .*
            sin.(last_phases[i] .+ 2pi * (0:buf_size-1) * freqs[i] / sample_rate)
        buf[:, 2] = buf[:, 2] .+ amp_mapping(u[i, 2]) .*
            sin.(last_phases[i] .+ 2pi * (0:buf_size-1) * freqs[i] / sample_rate)
        last_phases[i] = (last_phases[i] + 2pi * buf_size * freqs[i] / sample_rate) % 2pi
    end
    buf, last_phases
end

Base.@kwdef mutable struct AudioGen
        sample_rate::Float64
        fps::Int64
        freqs::Vector{Float64}
        last_phases = zeros(Float64, length(freqs))
        grain_time::Float64 = Float64(1 / fps)
        buf_size::Int64 = ceil(Int64, grain_time * sample_rate)
        buf::Matrix{Float64} = zeros(Float64, buf_size, 2)
end

# function setSineStack!(ag::AudioGen, buf::Vector{Float64}, u::Vector{Float64}, amp_mod::Function)
#     # mono version
#     for i = 1:length(freqs)
#         @inbounds @fastmath buf = buf .+ amp_mod(u[i]) .*
#             sin.(ag.last_phases[i] .+ 2pi * (0:ag.buf_size-1) * ag.freqs[i] / ag.sample_rate)
#         @inbounds @fastmath ag.last_phases[i] = (ag.last_phases[i] + 2pi * ag.buf_size * ag.freqs[i] / ag.sample_rate) % 2pi
#     end
# end

function setSineStack!(ag::AudioGen, u::Matrix{Float64}, amp_mod::Function)
    # stereo version
    ag.buf = zeros(Float64, ag.buf_size, 2) # reset buffer
    for i = 1:length(ag.freqs)
        @inbounds @fastmath ag.buf[:, 1] = ag.buf[:, 1] .+ amp_mod(u[i, 1]) .*
            sin.(ag.last_phases[i] .+ 2pi * (0:ag.buf_size-1) * ag.freqs[i] / ag.sample_rate)
        @inbounds @fastmath ag.buf[:, 2] = ag.buf[:, 2] .+ amp_mod(u[i, 2]) .*
            sin.(ag.last_phases[i] .+ 2pi * (0:ag.buf_size-1) * ag.freqs[i] / ag.sample_rate)
        @inbounds @fastmath ag.last_phases[i] = (ag.last_phases[i] + 2pi * ag.buf_size * ag.freqs[i] / ag.sample_rate) % 2pi
    end
end

function start(ks::KSIntegrator, audiogen::AudioGen,
        Nt::Int64, nplot::Int64, att::Float64, amp_mod::Function)
    PortAudio.PortAudioStream(0, 2; samplerate=sample_rate) do stream
        # initialize KS data and parameters
        # nsave = round(Int64, Nt/nplot)+1 # total number of presently displayed time steps
        # t = ks.dt*nplot * collect(0:nsave-1)
        # x = ks.Lx/ks.Nx * collect(1:ks.Nx)
        # U = zeros(Float64, nsave, ks.Nx)
        # U[end, :] = ks.u0

        # initialize figure and data
        U, t, x = ksIntegrate!(ks, Nt, nplot)
        U_max = maximum(abs.(U)) # for attenuation

        fig = GLMakie.Figure()
        ax = GLMakie.Axis(fig[1, 1])
        GLMakie.hidedecorations!(ax)
        hm = GLMakie.heatmap!(t, x, U,
            fxaa = true, inspectable = false,
            colorrange = (-U_max, U_max))
        GLMakie.display(fig)

        # live time stepping and plotting
        while GLMakie.events(fig).window_open[]
            evolve!(ks, nplot)
            U_new = getSolution(ks)
            U = vcat(U[begin+1:end,:], U_new')
            hm[3] = U

            # stereorize and normalize
            U_stereo = hcat(U_new, zeros(Float64, Nx))
            for i in 1:Nx
                if U_stereo[i] > 0
                    U_stereo[i, 2] = U_stereo[i, 1]
                    U_stereo[i, 1] = 0.0
                end
            end
            U_stereo = U_stereo .* att ./ U_max

            setSineStack!(audiogen, U_stereo, amp_mod)
            write(stream, audiogen.buf)
        end
    end
end


### main ###
# KS parameters
Lx = 128
Nx = 256
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
freqs = [40.0 + tanh(0.5*i/Nx) * 8000.0 for i in 0:Nx-1]
audiogen = AudioGen(sample_rate = sample_rate, fps = fps, freqs = freqs)

# start time stepping
Nt = 800
nplot = 1
att = 0.6 # additional attenuation
amp_mod = x -> x^5
start(ks, audiogen, Nt, nplot, att, amp_mod)
