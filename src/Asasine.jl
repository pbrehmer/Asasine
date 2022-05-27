module Asasine

using FFTW, GLMakie, PortAudio

include("./audio.jl")
export AudioGen, sine_stack!, stereorize

include("./ks.jl")
export KSIntegrator, evolve!, get_solution, integrate

include("./stepping.jl")
export start_stream 

end
