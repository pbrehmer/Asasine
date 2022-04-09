
Base.@kwdef mutable struct KSIntegrator
    u0::Vector{Float64}
    Lx::Float64
    dt::Float64

    ### additional physical parameters
    Nx::Integer = length(u0)
    kx::Vector{Integer} = vcat(0:Nx/2-1, 0, -Nx/2+1:-1) # shifted PBC integer wavenumbers
    alpha::Vector{Float64} = 2pi/Lx * kx

    # convenience variables
    L::Vector{Float64} = alpha.^2 - alpha.^4
    G::Vector{ComplexF64} = -0.5im * alpha
    dt2::Float64 = dt/2
    dt32::Float64 = 3*dt/2
    A::Vector{Float64} = ones(Nx) + dt2*L
    B::Vector{Float64} = (ones(Nx) - dt2*L).^-1
    
    # compute in-place planned FFTs
    FFT! = FFTW.plan_fft!((1+0im)*u0, flags=FFTW.ESTIMATE)
    IFFT! = FFTW.plan_ifft!((1+0im)*u0, flags=FFTW.ESTIMATE)

    # prepare time evolution in Fourier domain
    Nn::Vector{ComplexF64} = G .* FFTW.fft(u0.^2)
    Nn1::Vector{ComplexF64} = copy(Nn)
    u::Vector{ComplexF64} = FFT! * u0
end

function evolve!(ks::KSIntegrator, steps::Int=1)
    for n = 0:steps-1
        # compute Nn = G .* fft(real(ifft(u)).^2) inbounds with fastmath
        for i = 1:length(ks.Nn)
            @inbounds ks.Nn1[i] = ks.Nn[i]
            @inbounds ks.Nn[i] = ks.u[i]
        end

        ks.IFFT! * ks.Nn

        for i = 1:length(ks.Nn)
            @fastmath @inbounds ks.Nn[i] = ks.Nn[i] * ks.Nn[i]
        end
        
        ks.FFT! * ks.Nn
        
        for i = 1:length(ks.Nn)
            @fastmath @inbounds ks.Nn[i] = ks.G[i] * ks.Nn[i]
        end

        # compute u = B .* (A .* u + dt32*Nn - dt2*Nn1) inbounds with fastmath
        for i = 1:length(ks.u)
            @fastmath @inbounds ks.u[i] = ks.B[i] *
                (ks.A[i] * ks.u[i] + ks.dt32*ks.Nn[i] - ks.dt2*ks.Nn1[i])
        end
    end
end

function get_solution(ks::KSIntegrator)
    return real(FFTW.ifft(ks.u))
end

function integrate(ks::KSIntegrator, Nt::Integer, nplot::Integer)
    nsave = round(Integer, Nt/nplot)+1 # total number of saved time steps
    x = ks.Lx/ks.Nx * collect(1:ks.Nx)
    t = ks.dt*nplot * collect(0:nsave-1)
    U = zeros(Float64, nsave, ks.Nx)
    U[1,:] = ks.u0

    evolve!(ks)
    for i in 2:nsave
        U[i, :] = get_solution(ks)
        evolve!(ks, nplot)
    end

    return U, t, x
end
