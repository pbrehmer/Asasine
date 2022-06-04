
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

        # KS heatmap
        fig = Figure(; resolution=(1000, 700))
        ax1 = Axis(fig[1:3, 1]; title="Kuramoto-Sivashinsky stream")
        hidedecorations!(ax1)
        hm = heatmap!(
            ax1,
            t,
            x,
            U;
            colormap=:hawaii,
            fxaa=true,
            inspectable=false,
            colorrange=(-U_max, U_max),
        )
        Colorbar(fig[4, 1], hm; ticks=([-U_max, 0.0, U_max], ["L", " ", "R"]), vertical=false)

        # frequency spectrum
        ax2 = Axis(
            fig[1, 2:3];
            title="spectrum analyzer",
            xlabel="frequency [Hz]",
            xscale=log10,
            xticks=(
                [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000],
                string.([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]),
            ),
            ylabel="amplitude",
            yaxisposition=:right,
        )
        ylims!(ax2, 0.0, 0.75)
        xlims!(ax2, 20, audiogen.freqs[end] > 10000 ? audiogen.freqs[end] : 10000)
        hideydecorations!(ax2; label=false)
        spectrum_left = Observable{Vector{Float64}}(abs.(U_stereo[audiogen.freq_idx, 1]))
        spectrum_right = Observable{Vector{Float64}}(U_stereo[audiogen.freq_idx, 2])
        rb_left = rangebars!(
            ax2,
            audiogen.freqs,
            zeros(length(audiogen.freqs)),
            spectrum_left;
            color=:deeppink,
            label="L",
        )
        rb_right = rangebars!(
            ax2,
            audiogen.freqs,
            zeros(length(audiogen.freqs)),
            spectrum_right;
            color=:turquoise,
            label="R",
        )
        axislegend(ax2)

        # interactive elements
        sg = SliderGrid(
            fig[2, 2:3],
            (label="volume", range=0:0.1:1, startvalue=att),
            (label="fps", range=5:5:60, startvalue=audiogen.fps),
        )

        tb = Textbox(fig[3, 3]; placeholder="...", tellwidth=false, halign=:left)
        Label(fig[3, 2], "freq(x) = ")
        tbbutton = Button(fig[3, 3]; label="commit", tellwidth=false, halign=:right)

        running_toggle = Toggle(
            fig[4, 2]; active=true, height=20, tellwidth=false, halign=:left
        )
        Label(
            fig[4, 2],
            lift(x -> x ? "running" : "stopped", running_toggle.active);
            tellwidth=false,
            halign=:right,
            width=10.0,
        )

        # resize and display figure
        colsize!(fig.layout, 1, Relative(3 / 5))
        display(fig)

        # live time stepping and plotting
        while events(fig).window_open[]
            # write audio buffer to stream
            push = @task begin
                write(stream, audiogen.buf)
            end
            schedule(push)

            # generate new KS data
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
            schedule(generate)

            # update spectrum
            spectrum_left[] = abs.(U_stereo[audiogen.freq_idx, 1])
            spectrum_right[] = U_stereo[audiogen.freq_idx, 2]

            # get interactive updates
            on(sg.sliders[1].value) do volume
                att = volume
            end
            on(sg.sliders[2].value) do fps
                audiogen = AudioGen(
                    audiogen.sample_rate, fps, audiogen.freq_func; freq_idx=audiogen.freq_idx
                )
            end
            # on(tb.stored_string) do s
            #     func_string[] = s
            # end
            # on(func_string) do val
            #     eval(Meta.parse("f(x) = " * val))
            #     audiogen.freq_func = f
            #     audiogen.freqs = audiogen.freq_func.(Float64.(audiogen.freq_idx))
            # end
            on(tbbutton.clicks) do n
                # eval(Meta.parse("parsedfunc(x) = " * tb.stored_string.val))
                # audiogen = AudioGen(
                #     audiogen.sample_rate, audiogen.fps, parsedfunc; freq_idx=audiogen.freq_idx
                # )
            end

            wait(push)
            while !running_toggle.active.val && events(fig).window_open[]
                sleep(0.2)
            end
        end
    end
end
