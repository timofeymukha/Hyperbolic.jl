module FluxM

export Flux, FluxBurgers, FluxLinear, FluxBuckleyLeverett, FluxTraffic,
       FluxTrafficMultilane,
       flux, flux_deriv

abstract Flux

type FluxBurgers <: Flux
    qS::Real
    function FluxBurgers()
        qS = 0
        new(qS)
    end
end

type FluxBuckleyLeverett <: Flux
    qS::Real
    function FluxBuckleyLeverett()
        qS = 1
        new(qS)
    end
end

type FluxLinear <: Flux
    a::Real
end

type FluxTraffic <:Flux
    uMax::Real
    qS::Real

    function FluxTraffic(uMax)
        qS = 0.5
        new(uMax, qS)
    end
end

type FluxTrafficMultilane <:Flux
    uMax::Real
    qS::Real
    a::Array

    function FluxTrafficMultilane(uMax, a)
        qS = 0.5
        new(uMax, qS, a)
    end
end

# Traffic flow with variable lanes
function flux(fluxType::FluxTrafficMultilane, q)
    return fluxType.uMax*q.*(1 - q./fluxType.a)
end

function flux_deriv(fluxType::FluxTrafficMultilane, q)
    a = fluxType.a
    return fluxType.uMax.*(1 - 2*q./a -q.^2./a.^2)
end

# Traffic flow flux
function flux(fluxType::FluxTraffic, q)
    return fluxType.uMax*q.*(1 - q)
end

function flux_deriv(fluxType::FluxTraffic, q)
    return fluxType.uMax.*(1 - 2*q)
end

# Burgers equation flux
function flux(fluxType::FluxBurgers, q)
    return 0.5*q.^2
end

function flux_deriv(fluxType::FluxBurgers, q)
    return 1.0*q
end

# Simple linear advection flux
function flux(fluxType::FluxLinear, q)
    return fluxType.a*q
end

function flux_deriv(fluxType::FluxLinear, q)
    return fluxType.a*ones(Float64, size(q))
end

# Byckley-Leverett flux
function flux(fluxType::FluxBuckleyLeverett, q)
    return q.^2./(q.^2 + 0.5*(1 - q).^2)
end

function flux_deriv(fluxType::FluxBuckleyLeverett, q)
    return -4*(q-1).*q./(3*q.^2 - 2*q + 1).^2
end



end
