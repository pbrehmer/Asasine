# Open an audio stream on the system's default output and
# live animate and sonify the KS solution while time stepping
function stream(sbuf::SampleBuf{T}, son::Sonifier{T}, ks::KSIntegrator;
                Nt=600, nplot=1, att=0.8, resolution=(1000, 600), theme=kstheme) where {T}
    PortAudioStream(0, 2; samplerate=samplerate(sbuf)) do stream
        # GLMakie settings
        GLMakie.activate!(inline=false)  # Open visuals in window
        set_theme!(theme)

        # Initialize frame
        U, t, x = integrate(ks, Nt; nplot)
        Umax = maximum(abs, U)  # For attenuation
        Ustereo = stereorize(U[end, :])

        # KS heatmap
        fig = Figure(; resolution)
        ax1 = Axis(fig[1:4, 1]; title="Kuramoto-Sivashinsky stream")
        hidedecorations!(ax1)
        hm = heatmap!(ax1, t, x, U;
                      colormap=:hawaii, inspectable=false, colorrange=(-Umax, Umax))

        # Spectrum analyzer
        ax2 = Axis(fig[1, 2];
                   title="spectrum analyzer", xlabel="frequency [Hz]", xscale=log10,
                   xticks=([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000],
                           string.([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]),),
                   ylabel="amplitude", yaxisposition=:right)
        ylims!(ax2, 0.0, 0.8)
        xlims!(ax2, 20, maximum(son.freqs) > 10000 ? maximum(son.freqs) : 10000)
        hideydecorations!(ax2; label=false)
        spectruml = Observable{Vector{T}}(abs.(Ustereo[son.freqidx, 1]))
        spectrumr = Observable{Vector{T}}(abs.(Ustereo[son.freqidx, 2]))
        rangebars!(ax2, son.freqs, zeros(length(son.freqs)), spectruml;
                   color=:deeppink, label="L")
        rangebars!(ax2, son.freqs, zeros(length(son.freqs)), spectrumr;
                   color=:turquoise, label="R")
        axislegend(ax2)

        # UI elements
        sg = SliderGrid(fig[2, 2],
                        (; label="volume", range=0:0.1:1, startvalue=att),
                        (; label="fps", range=5:5:60, startvalue=1/domain(sbuf)[end]))
        running_toggle = Toggle(fig[3, 2];
                                active=true, height=20, tellwidth=false, halign=0.225)
        Label(fig[3, 2], lift(x -> x ? "running" : "stopped", running_toggle.active);
              tellwidth=false, halign=:left, width=10.0)
        Colorbar(fig[4, 2], hm;
                 ticks=([-Umax, 0.0, Umax], ["L", " ", "R"]), vertical=false)

        # Resize and display figure
        colsize!(fig.layout, 1, Relative(3 / 5))
        display(fig)

        # Live time stepping and plotting
        sbuf .= 0.0
        lastphases = zeros(T, length(son.freqs))
        while events(fig).window_open[]
            # Write audio buffer to stream
            write(stream, sbuf)

            # Generate new KS data
            evolve!(ks, nplot)
            Unew = solution(ks)
            U = vcat(U[2:end, :], Unew')
            hm[3] = U

            # Stereorize and normalize
            Ustereo .= stereorize(Unew)
            Unewmax = maximum(abs.(Unew))
            if Unewmax > Umax  # Update normalization if necessary
                Umax = Unewmax
            end
            Ustereo = Ustereo .* T(att) ./ Umax
            sinestack!(sbuf, Ustereo, son, lastphases)

            # Update spectrum
            spectruml[] = abs.(Ustereo[son.freqidx, 1])
            spectrumr[] = abs.(Ustereo[son.freqidx, 2])

            # Get interactive updates
            on(sg.sliders[1].value) do volume
                att = T(volume)
            end
            on(sg.sliders[2].value) do fps
                sbuf = SampleBuf(T, samplerate(sbuf), (1 / fps)s, nchannels(sbuf))
            end

            # Halt live stepping
            while !running_toggle.active.val && events(fig).window_open[]
                sleep(0.2)
            end
        end
    end
end
