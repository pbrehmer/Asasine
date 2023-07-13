module Asasine

using FFTW
using GLMakie
using PortAudio
using SampledSignals
using LibSndFile

include("theme.jl")
export kstheme

include("./audio.jl")
export Sonifier, sinestack!, stereorize

include("./ks.jl")
export KSIntegrator, evolve!, solution, integrate

include("./stream.jl")
export stream

include("render.jl")
export renderaudio, renderimage, rendervideo

end
