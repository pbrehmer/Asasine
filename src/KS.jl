
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

function get_solution(ks::KSIntegrator)
    return real(FFTW.ifft(ks.u))
end

function integrate(ks::KSIntegrator, Nt::Integer, nplot::Integer)
    nsave = round(Integer, Nt/nplot)+1 # total number of saved time steps
    x = ks.Lx/ks.Nx * collect(1:ks.Nx)
    t = ks.dt*nplot * collect(0:nsave-1)
    U = zeros(Float64, nsave, ks.Nx)
    U[1,:] = ks.u0

    evolve!(ks)
    for i in 2:nsave
        U[i, :] = get_solution(ks)
        evolve!(ks, nplot)
    end

    return U, t, x
end

function start(ks::KSIntegrator, audiogen::AudioGen,
        Nt::Integer, nplot::Integer, att::Float64, amp_mod::Function)
    PortAudio.PortAudioStream(0, 2; samplerate=sample_rate) do stream
        # initialize KS data and parameters
        # nsave = round(Integer, Nt/nplot)+1 # total number of presently displayed time steps
        # t = ks.dt*nplot * collect(0:nsave-1)
        # x = ks.Lx/ks.Nx * collect(1:ks.Nx)
        # U = zeros(Float64, nsave, ks.Nx)
        # U[end, :] = ks.u0

        # initialize figure and data
        U, t, x = integrate(ks, Nt, nplot)
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
            U_new = get_solution(ks)
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

            sine_stack!(audiogen, U_stereo, amp_mod)
            write(stream, audiogen.buf)
        end
    end
end