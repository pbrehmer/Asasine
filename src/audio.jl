
mutable struct AudioGen
    sample_rate::Float64
    fps::Int
    freq_func::Function
    freq_idx::StepRange
    freqs::Vector{Float64}
    last_phases::Vector{Float64}
    grain_time::Float64
    buf_size::Int
    buf::Matrix{Float64}
end
function AudioGen(sample_rate, fps, freq_func::Function; freq_idx=1:1:length(freqs))
    freqs = freq_func.(Float64.(freq_idx))
    last_phases = zeros(Float64, length(freqs))
    grain_time::Float64 = Float64(1 / fps)
    buf_size::Integer = ceil(Integer, grain_time * sample_rate)
    buf::Matrix{Float64} = zeros(Float64, buf_size, 2)

    AudioGen(sample_rate, fps, freq_func, freq_idx, freqs, last_phases, grain_time, buf_size, buf)
end

function sine_stack!(ag::AudioGen, u::Matrix{Float64}, amp_mod::Function)
    # stereo version
    ag.buf = zeros(Float64, ag.buf_size, 2) # reset buffer
    for i in 1:length(ag.freqs)
        @inbounds @fastmath ag.buf[:, 1] = ag.buf[:, 1] .+ amp_mod(u[ag.freq_idx[i], 1]) .*
            sin.(ag.last_phases[i] .+ 2pi * (0:ag.buf_size-1) * ag.freqs[i] / ag.sample_rate)
        @inbounds @fastmath ag.buf[:, 2] = ag.buf[:, 2] .+ amp_mod(u[ag.freq_idx[i], 2]) .*
            sin.(ag.last_phases[i] .+ 2pi * (0:ag.buf_size-1) * ag.freqs[i] / ag.sample_rate)
        @inbounds @fastmath ag.last_phases[i] = (ag.last_phases[i] + 2pi * ag.buf_size * ag.freqs[i] / ag.sample_rate) % 2pi
    end
end

function stereorize(buffer_mono::Vector{Float64})
    buffer_stereo = hcat(buffer_mono, zeros(Float64, length(buffer_mono)))
    for i = 1:length(buffer_mono)
        @inbounds if buffer_stereo[i] > 0
            @inbounds buffer_stereo[i, 2] = buffer_stereo[i, 1]
            @inbounds buffer_stereo[i, 1] = 0.0
        end
    end

    buffer_stereo
end
