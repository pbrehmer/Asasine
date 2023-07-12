# Contains settings and convenience variables that are used to integrate the KS equation
# T: floating-point type (for audio operations usually 32 bit floats)
# C: complex floating point type
@kwdef mutable struct KSIntegrator{T<:Real,C}
    u0::Vector{T}  # Initial condition
    Lx::T  # Spatial length
    dt::T  # Time step size

    Nx::Int = length(u0)  # Spatial resolution
    kx::Vector{T} = vcat(0:Nx/2-1, 0, -Nx/2+1:-1)  # PBC wave numbers
    alpha::Vector{T} = 2pi / Lx * kx  # Converted wave numbers

    # Convenience variables
    L::Vector{T} = alpha.^2 - alpha.^4
    G::Vector{C} = -0.5im * alpha
    dt2::T = dt / 2
    dt32::T = 3 * dt / 2
    A::Vector{T} = ones(Nx) + dt2 * L
    B::Vector{T} = (ones(Nx) - dt2 * L).^-1

    # Compute in-place planned FFTs
    FFT! = plan_fft!((1 + 0im) * u0)
    IFFT! = plan_ifft!((1 + 0im) * u0)

    # Time evolution in Fourier domain
    Nn::Vector{C} = G .* fft(u0 .^ 2)
    Nn1::Vector{C} = deepcopy(Nn)
    u::Vector{C} = fft(u0)
end

# Default to 32 bit floating point arithmetic
KSIntegrator{T}(u0, Lx, dt) where {T<:Real} = KSIntegrator{T,Complex{T}}(; u0=T.(u0), Lx=T(Lx), dt=T(dt))

# Obtain real-space solution
solution(ks::KSIntegrator) = real(ifft(ks.u))

# Perform one CNAB2 time evolution step
function evolve!(ks::KSIntegrator, steps::Int=1)
    for _ = 0:steps-1
        # Compute Nn = G .* fft(real(ifft(u)).^2) inbounds with fastmath
        for i in 1:length(ks.Nn)
            @inbounds ks.Nn1[i] = ks.Nn[i]
            @inbounds ks.Nn[i] = ks.u[i]
        end
        ks.IFFT! * ks.Nn
        for i in 1:length(ks.Nn)
            @fastmath @inbounds ks.Nn[i] = ks.Nn[i] * ks.Nn[i]
        end
        ks.FFT! * ks.Nn
        for i in 1:length(ks.Nn)
            @fastmath @inbounds ks.Nn[i] = ks.G[i] * ks.Nn[i]
        end

        # Compute u = B .* (A .* u + dt32*Nn - dt2*Nn1) inbounds with fastmath
        for i in 1:length(ks.u)
            @fastmath @inbounds ks.u[i] =
                ks.B[i] * (ks.A[i] * ks.u[i] + ks.dt32 * ks.Nn[i] - ks.dt2 * ks.Nn1[i])
        end
    end
end

# Integrate KS equation for Nt time steps, saving the solution each nplot steps
function integrate(ks::KSIntegrator{T}, Nt; nplot=1) where {T}
    nsave = round(Integer, Nt / nplot) + 1  # Total number of saved time steps
    x = ks.Lx / ks.Nx * collect(1:ks.Nx)
    t = ks.dt * nplot * collect(0:nsave-1)
    U = zeros(T, nsave, ks.Nx)
    U[1, :] = ks.u0

    for i in 1:nsave
        evolve!(ks, nplot)
        U[i, :] = solution(ks)
    end

    U, t, x
end
