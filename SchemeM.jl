module SchemeM

using Hyperbolic.FluxM
using Hyperbolic.LimiterM
using Hyperbolic.MeshM

export Scheme, LaxFriedrichs, Godunov, Roe, UpwindSimple, LocalLaxFriedrichs,
numerical_flux, flux_correction, accepts_limiter

abstract Scheme

type LaxFriedrichs <: Scheme
    flux::Flux
end

type LocalLaxFriedrichs <: Scheme
    flux::Flux
end

type UpwindSimple <: Scheme
    flux::Flux
end

type Godunov <: Scheme
    flux::Flux
    limiter::Limiter
end

type Roe <: Scheme
    flux::Flux
end

function numerical_flux(mesh::Mesh1d, scheme::UpwindSimple, Q)
    APlus = zeros(Float64, length(Q)+1)
    AMinus = zeros(Float64, length(Q)+1)

    for i in 2:size(Q, 1)-1
        if Q[i] >= 0
            APlus[i] = Q[i]*Q[i]
            AMinus[i+1] = -Q[i-1]*Q[i]
        else

            APlus[i] = Q[i]*Q[i+1]
            AMinus[i+1] = -Q[i]*Q[i]
        end
    end

    return APlus, AMinus, zeros(Float64, size(APlus))
end

function numerical_flux(mesh::Mesh1d, scheme::LaxFriedrichs, Q)
    F = zeros(Float64, length(Q)+1)

    fluxType = scheme.flux
    fluxValue = Vector{Float64}(flux(fluxType, Q))

    a = mesh.h/mesh.dt

    F[2:end-1] = 0.5*(fluxValue[1:end-1] +fluxValue[2:end] -
                      a*(Q[2:end] - Q[1:end-1]))

    F[1] = F[2]
    F[end] = F[end-1]

    return -F, F, zeros(Float64, size(F))
end

function numerical_flux(mesh::Mesh1d, scheme::LocalLaxFriedrichs, Q)
    F = zeros(Float64, length(Q)+1)
    a = zeros(Float64, length(Q)+1)

    fluxType = scheme.flux
    fluxValue = Vector{Float64}(flux(fluxType, Q))
    fluxDerivValue = Vector{Float64}(flux_deriv(fluxType, Q))

    for i in 2:size(a, 1)-1
        a[i] = maximum([abs(fluxDerivValue[i-1]), abs(fluxDerivValue[i])])
    end

    F[2:end-1] = 0.5*(fluxValue[1:end-1] +fluxValue[2:end] -
                      a[2:end-1].*(Q[2:end] - Q[1:end-1]))

    F[1] = F[2]
    F[end] = F[end-1]

    return -F, F, zeros(Float64, size(F))
end

function numerical_flux(mesh::Mesh1d, scheme::Roe, Q)
    F = zeros(Float64, length(Q)+1)
    a = zeros(Float64, length(Q)+1)

    fluxValue = Vector{Float64}(flux(scheme.flux, Q))

    for i in 2:size(a, 1)-1
        if Q[i] != Q[i-1]
            a[i] = abs((fluxValue[i] - fluxValue[i-1])/(Q[i]- Q[i-1]))
        end
    end

    F[2:end-1] = 0.5*(fluxValue[1:end-1] +fluxValue[2:end] -
                      a[2:end-1].*(Q[2:end] - Q[1:end-1]))

    F[1] = F[2]
    F[end] = F[end-1]

    return -F, F, zeros(Float64, size(F))
end

function numerical_flux(mesh::Mesh1d, scheme::Godunov, Q)

    W = zeros(Float64, length(Q)+1)
    s = zeros(Float64, length(Q)+1)
    APlus = zeros(Float64, length(Q)+1)
    AMinus = zeros(Float64, length(Q)+1)
    FTilde = zeros(Float64, length(Q)+1)

    fluxType = scheme.flux
    fluxValue = Vector{Float64}(flux(fluxType, Q))
    fluxDerivValue = Vector{Float64}(flux_deriv(fluxType, Q))

    # Nemann conidtions
    insert!(fluxValue, 1, fluxValue[1])
    insert!(fluxDerivValue, 1, fluxDerivValue[1])
    push!(fluxValue, fluxValue[end])
    push!(fluxDerivValue, fluxDerivValue[end])

    W[2:end-1] = Q[2:end]-Q[1:end-1]

    # Neumann boundaries
    W[1] = 0
    W[end] = 0

    for i in eachindex(s)
        if W[i] != 0
            s[i] = (fluxValue[i+1] - fluxValue[i])/W[i]
        end
    end


    # See (12.8) in LeVeque
    for i in eachindex(APlus)
        APlus[i] = max(s[i], 0)*W[i]
        AMinus[i] = min(s[i], 0)*W[i]
    end

    # Entropy fix, see (12.9) in LeVeque
    for i in eachindex(APlus)
        if fluxDerivValue[i] < 0 && fluxDerivValue[i+1] > 0
    #        APlus[i] = fluxValue[i+1] - flux(fluxType, fluxType.qS)
    #        AMinus[i] = flux(fluxType, fluxType.qS) - fluxValue[i]
        end
    end

    FTilde = flux_correction(mesh, scheme.limiter, W, s)
    return APlus, AMinus, FTilde
end

function flux_correction(mesh::Mesh1d,
                         limiterType::Limiter,
                         W::Vector{Float64},
                         s::Vector{Float64})
    theta = zeros(size(W))

    for i in 2:size(theta, 1)-1
        if s[i] >= 0 && W[i] != 0
            theta[i] = W[i-1]/W[i]
        elseif W[i] != 0
            theta[i] = W[i+1]/W[i]
        end
    end

    theta[1] = theta[2]

    for i in eachindex(W)
        W[i] = W[i] * limiter(limiterType, theta[i])
    end


    F = 0.5*abs(s).*(1-mesh.dt/mesh.h*abs(s)).*W

    return F
end

end
