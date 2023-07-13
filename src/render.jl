function renderaudio(sbuf::SampleBuf{T}, son::Sonifier{T}, ks::KSIntegrator,
                     filename::AbstractString; Nt=600, nplot=1) where {T}
    U, = integrate(ks, Nt; nplot)

    sbuftot = deepcopy(sbuf)  # To not mutate sbuf
    datatot = Matrix{eltype(sbuf)}(undef, 0, nchannels(sbuf))
    lastphases = zeros(T, length(son.freqs))
    for u in eachslice(U; dims=1)
        ustereo = stereorize(u)
        sinestack!(sbuftot, ustereo, son, lastphases)
        datatot = vcat(datatot, sbuftot.data)
    end
    stereomax = maximum(sum(datatot; dims=2))
    sbuftot.data = datatot ./ stereomax  # Normalize volume to 1

    save(filename, sbuftot)
    sbuftot
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