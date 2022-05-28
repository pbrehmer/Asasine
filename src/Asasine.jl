module Asasine

using FFTW, GLMakie, PortAudio

black = :gray10
green = :seagreen1
green_dimmed = :seagreen
green_lessdimmed = :mediumseagreen
colormap = :hawaii
ks_theme = Theme(;
    backgroundcolor=black,
    textcolor=green,
    linecolor=:white,
    font="FreeMono Bold",
    fontsize=14,
    Axis=(
        titlefont="FreeMono Bold",
        titlesize=16,
        backgroundcolor=:transparent,
        bottomspinecolor=green,
        topspinecolor=green,
        leftspinecolor=green,
        rightspinecolor=green,
        xgridcolor=RGBAf(1, 1, 1, 0.16),
        ygridcolor=RGBAf(1, 1, 1, 0.16),
        xtickcolor=green,
        ytickcolor=green,
    ),
    Legend=(framecolor=green, bgcolor=black),
    Slider=(
        color_active=green,
        color_active_dimmed=green_dimmed,
        color_inactive=green_dimmed,
    ),
    Colorbar=(
        tickcolor=green,
        spinecolor=green,
        topspinecolor=green,
        bottomspinecolor=green,
        leftspinecolor=green,
        rightspinecolor=green,
        colormap=colormap,
    ),
    Toggle = (
        buttoncolor=green,
        toggleduration=0.1,
        # framecolor_active=:white,
        # framecolor_inactive=green_dimmed,
    ),
)
export ks_theme

include("./audio.jl")
export AudioGen, sine_stack!, stereorize

include("./ks.jl")
export KSIntegrator, evolve!, get_solution, integrate

include("./stepping.jl")
export start_stream

end
