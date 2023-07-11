
struct Sonifier{T<:Real}
    freqs::Vector{T}
    freqidx::AbstractVector{Int}
    ampmap
end

# Stereo
function sinestack!(sbuf::SampleBuf{T,2}, u::Matrix{T}, son::Sonifier, lastphases) where {T}
    sbuf .= 0.0  # Reset buffer
    dom = domain(sbuf)
    for (i, ν) in enumerate(son.freqs)
        @inbounds @fastmath sbuf[:, 1] .+= son.ampmap(u[son.freqidx[i], 1]) .*
            sin.(lastphases[i] .+ 2pi * dom * ν)
        @inbounds @fastmath sbuf[:, 2] .+= son.ampmap(u[son.freqidx[i], 2]) .*
            sin.(lastphases[i] .+ 2pi * dom * ν)
        @inbounds @fastmath lastphases[i] = (lastphases[i] + 2pi * dom[end] * ν ) % 2pi
    end
end

# Mono
function sinestack!(sbuf::SampleBuf{T,1}, u::Vector{T}, freqs, lastphases, ampmap) where {T}
    sbuf .= 0.0  # Reset buffer
    dom = domain(sbuf)
    for (i, ν) in enumerate(freqs)
        @inbounds @fastmath sbuf .+= ampmap(u[ag.freq_idx[i], 1]) .*
            sin.(lastphases[i] .+ 2pi * dom * ν)
        @inbounds @fastmath lastphases[i] = (lastphases[i] + 2pi * dom * ν ) % 2pi
    end
end

# TODO: sine stack generation via inverse FFT

function stereorize(u::Vector{T}) where {T}
    ustereo = hcat(u, zeros(T, length(u)))
    for i in 1:length(u)
        @inbounds if u[i] > 0
            @inbounds ustereo[i, 2] = ustereo[i, 1]
            @inbounds ustereo[i, 1] = 0.0
        end
    end

    ustereo
end
