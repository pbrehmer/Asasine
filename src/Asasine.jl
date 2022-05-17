module Asasine

using FFTW, GLMakie, PortAudio
import Makie: PriorityObservable
set_theme!(theme_black())

include("./audio.jl")
export AudioGen, sine_stack!

include("./ks.jl")
export KSIntegrator, evolve!, get_solution, integrate

include("./stepping.jl")
export start

end
