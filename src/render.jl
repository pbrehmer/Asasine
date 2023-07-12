function renderaudio(sbuf::SampleBuf{T}, son::Sonifier{T}, ks::KSIntegrator;
                     Nt=600, nplot=1, att=0.8) where {T}

end

function renderimage(ks::KSIntegrator, filename::AbstractString;
                     Nt=600, nplot=1, resolution=(1200, 800), theme=kstheme,
                     hmkwargs=(; colormap=:hawaii))
    set_theme!(theme)
    U, t, x = integrate(ks, Nt; nplot)
    Umax = maximum(abs, U)

    fig, ax, hm = heatmap(t, x, U; hmkwargs..., colorrange=(-Umax, Umax))
    hidedecorations!(ax)
    save(filename, fig; resolution)

    fig, ax, hm
end

function rendervideo(ks::KSIntegrator; Nt=600, nplot=1, resolution=(1000, 600), theme=kstheme)

end