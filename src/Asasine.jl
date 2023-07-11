module Asasine

using FFTW
using GLMakie
using PortAudio
using SampledSignals

include("theme.jl")
export kstheme

include("./audio.jl")
export Sonifier, sinestack!, stereorize

include("./ks.jl")
export KSIntegrator, evolve!, solution, integrate

include("./stepping.jl")
export stream

end
