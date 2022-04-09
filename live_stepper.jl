import FFTW, GLMakie, PortAudio
using Makie: PriorityObservable
GLMakie.set_theme!(GLMakie.theme_black())

include("./src/Audio.jl")
include("./src/KS.jl")

function start(ks::KSIntegrator, audiogen::AudioGen,
        Nt::Integer, nplot::Integer, att::Float64, amp_mod::Function)
    PortAudio.PortAudioStream(0, 2; samplerate=sample_rate) do stream
        # initialize KS data and parameters
        #=
        nsave = round(Integer, Nt/nplot)+1 # total number of presently displayed time steps
        t = ks.dt*nplot .* collect(0:nsave-1)
        x = ks.Lx/ks.Nx .* collect(1:ks.Nx)
        U = zeros(Float64, nsave, ks.Nx)
        U[end, :] = ks.u0
        =#

        # initialize figure and data
        U, t, x = integrate(ks, Nt, nplot)
        U_max = maximum(abs.(U)) # for attenuation

        fig = GLMakie.Figure()
        ax = GLMakie.Axis(fig[1, 1])
        GLMakie.hidedecorations!(ax)
        hm = GLMakie.heatmap!(fig[1, 1], t, x, U,
            fxaa = true, inspectable = false,
            colormap = :PRGn_11, colorrange = (-U_max, U_max))
        GLMakie.Colorbar(fig[1, 2], hm, ticks = ([-U_max, U_max], ["L", "R"]))
        GLMakie.display(fig)

        # live time stepping and plotting
        while GLMakie.events(fig).window_open[]
            evolve!(ks, nplot)
            U_new = get_solution(ks)
            U = vcat(U[begin+1:end,:], U_new')
            hm[3] = U

            # stereorize and normalize
            U_stereo = hcat(U_new, zeros(Float64, Nx))
            for i in 1:Nx
                @inbounds if U_stereo[i] > 0
                    @inbounds U_stereo[i, 2] = U_stereo[i, 1]
                    @inbounds U_stereo[i, 1] = 0.0
                end
            end
            U_new_max = maximum(abs.(U_new)) # update normalization if necessary
            if U_new_max > U_max
                U_max = U_new_max
            end
            U_stereo = U_stereo .* att ./ U_max

            sine_stack!(audiogen, U_stereo, amp_mod)
            write(stream, audiogen.buf)
        end
    end
end

# KS parameters
Lx = 80
Nx = 1024
dt = 1/32
x = Lx/Nx * collect(1:Nx)
u0 = cos.(x) .+ 0.1*cos.(x/16) .* (1 .+ 2*sin.(x/16));
# u0 = 3.0 * exp.(-(x .- Nx/2).^2 ./ sqrt(Nx))
# u0 = 2.0 * (rand(length(x)) .- 0.5)
# u0 = 2.0 * rand(length(x))
ks = KSIntegrator(u0 = u0, Lx = Lx, dt = dt);

# audio parameters
sample_rate = Float64(44100)
fps = 30
downsampler = 8
println("Video resolution per time step: ", Nx)
println("Audio resolution per time step: ", Nx / downsampler)
# freqs = collect(range(40, 400, Nx))
freqs = [50.0 + tanh(0.3*i/Nx) * 2500.0 for i in 0:downsampler:Nx-1]
audiogen = AudioGen(sample_rate = sample_rate, fps = fps, freqs = freqs);

# start time stepping
Nt = 800
nplot = 1
att = 0.6 # additional attenuation
amp_mod = x -> x^5
start(ks, audiogen, Nt, nplot, att, amp_mod)
