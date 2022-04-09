
#=
function sine_stack!(ag::AudioGen, buf::Vector{Float64}, u::Vector{Float64}, amp_mod::Function)
    # mono version
    for i = 1:length(freqs)
        @inbounds @fastmath buf = buf .+ amp_mod(u[i]) .*
            sin.(ag.last_phases[i] .+ 2pi * (0:ag.buf_size-1) * ag.freqs[i] / ag.sample_rate)
        @inbounds @fastmath ag.last_phases[i] = (ag.last_phases[i] + 2pi * ag.buf_size * ag.freqs[i] / ag.sample_rate) % 2pi
    end
end
=#

function sine_stack!(ag::AudioGen, u::Matrix{Float64}, amp_mod::Function)
    # stereo version
    ag.buf = zeros(Float64, ag.buf_size, 2) # reset buffer
    for i = 1:length(ag.freqs)
        @inbounds @fastmath ag.buf[:, 1] = ag.buf[:, 1] .+ amp_mod(u[i, 1]) .*
            sin.(ag.last_phases[i] .+ 2pi * (0:ag.buf_size-1) * ag.freqs[i] / ag.sample_rate)
        @inbounds @fastmath ag.buf[:, 2] = ag.buf[:, 2] .+ amp_mod(u[i, 2]) .*
            sin.(ag.last_phases[i] .+ 2pi * (0:ag.buf_size-1) * ag.freqs[i] / ag.sample_rate)
        @inbounds @fastmath ag.last_phases[i] = (ag.last_phases[i] + 2pi * ag.buf_size * ag.freqs[i] / ag.sample_rate) % 2pi
    end
end
