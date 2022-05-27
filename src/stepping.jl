
function start_stream(
    sample_rate, ks::KSIntegrator, audiogen::AudioGen, Nt, nplot, att, amp_mod::Function
)
    PortAudioStream(0, 2; samplerate=sample_rate) do stream
        # initialize KS data and parameters
        # nsave = round(Integer, Nt/nplot)+1 # total number of presently displayed time steps
        # t = ks.dt*nplot .* collect(0:nsave-1)
        # x = ks.Lx/ks.Nx .* collect(1:ks.Nx)
        # U = zeros(Float64, nsave, ks.Nx)
        # U[end, :] = ks.u0

        # initialize data
        U, t, x = integrate(ks, Nt, nplot)
        U_max = maximum(abs.(U)) # for attenuation
        U_stereo = stereorize(U[end, :])

        # set up heatmap
        fig = Figure(; resolution=(1000, 700))
        ax1 = Axis(fig[1:2, 1])
        hidedecorations!(ax1)
        hm = heatmap!(
            ax1,
            t,
            x,
            U;
            fxaa=true,
            inspectable=false,
            colormap=:hawaii,
            colorrange=(-U_max, U_max),
        )
        Colorbar(fig[3, 1], hm; ticks=([-U_max, U_max], ["L", "R"]), vertical=false)

        # set up frequency spectrum
        ax2 = Axis(fig[1, 2]; xscale=log10)
        ylims!(ax2, 0.0, 0.75)
        freq_step = ks.Nx ÷ length(audiogen.freqs)
        freq_idx = 1:freq_step:ks.Nx
        spectrum_left = Observable{Vector{Float64}}(abs.(U_stereo[freq_idx, 1]))
        spectrum_right = Observable{Vector{Float64}}(U_stereo[freq_idx, 2])

        rangebars!(
            ax2,
            audiogen.freqs,
            zeros(length(audiogen.freqs)),
            spectrum_left;
            color=:magenta,
        )
        rangebars!(
            ax2,
            audiogen.freqs,
            zeros(length(audiogen.freqs)),
            spectrum_right;
            color=:cyan,
        )

        display(fig)

        # live time stepping and plotting
        while events(fig).window_open[]
            push = @task begin
                write(stream, audiogen.buf)
            end
            schedule(push)
            generate = @task begin
                evolve!(ks, nplot)
                U_new = get_solution(ks)
                U = vcat(U[begin+1:end, :], U_new')
                hm[3] = U

                # stereorize and normalize
                U_stereo .= stereorize(U_new)
                U_new_max = maximum(abs.(U_new))
                if U_new_max > U_max # update normalization if necessary
                    U_max = U_new_max
                end
                U_stereo = U_stereo .* att ./ U_max
                sine_stack!(audiogen, U_stereo, amp_mod)
            end
            spectrum_left[] = abs.(U_stereo[freq_idx, 1])
            spectrum_right[] = U_stereo[freq_idx, 2]

            schedule(generate)
            wait(push)
        end
    end
end
