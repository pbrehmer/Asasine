
function start(
    sample_rate,
    ks::KSIntegrator,
    audiogen::AudioGen,
    Nt,
    nplot,
    att,
    amp_mod::Function,
)
    PortAudioStream(0, 2; samplerate=sample_rate) do stream
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

        fig = Figure()
        ax = Axis(fig[1, 1])
        hidedecorations!(ax)
        hm = heatmap!(
            fig[1, 1],
            t,
            x,
            U;
            fxaa=true,
            inspectable=false,
            colormap=:PRGn_11,
            colorrange=(-U_max, U_max),
        )
        Colorbar(fig[1, 2], hm; ticks=([-U_max, U_max], ["L", "R"]))
        display(fig)

        # live time stepping and plotting
        while events(fig).window_open[]
            evolve!(ks, nplot)
            U_new = get_solution(ks)
            U = vcat(U[begin+1:end, :], U_new')
            hm[3] = U

            # stereorize and normalize
            U_stereo = stereorize(U_new)
            U_new_max = maximum(abs.(U_new))
            if U_new_max > U_max # update normalization if necessary
                U_max = U_new_max
            end
            U_stereo = U_stereo .* att ./ U_max

            sine_stack!(audiogen, U_stereo, amp_mod)
            write(stream, audiogen.buf)
        end
    end
end
