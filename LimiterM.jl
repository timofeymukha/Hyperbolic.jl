module LimiterM

export Limiter, Upwind, LaxWendroff, BeamWarming, Fromm, VanLeer, MC,
Superbee, Minmod, limiter

abstract Limiter

type Upwind <: Limiter
end

type LaxWendroff <: Limiter
end

type BeamWarming <: Limiter
end

type Fromm <: Limiter
end

type VanLeer <: Limiter
end

type MC <: Limiter
end

type Superbee <: Limiter
end

type Minmod <: Limiter
end

function limiter(limiterType::Upwind, theta::Real)
    return 0.0
end

function limiter(limiterType::LaxWendroff, theta::Real)
    return 1.0
end

function limiter(limiterType::BeamWarming, theta::Real)
    return 1*theta
end

function limiter(limiterType::Fromm, theta::Real)
    return 0.5*(1 + theta)
end

function limiter(limiterType::VanLeer, theta::Real)
    return (theta + abs(theta))/(1+ abs(theta))
end



function limiter(limiterType::MC, theta::Real)
    return max(0.0, min(0.5*(1 + theta), 2.0, 2*theta))
end

function limiter(limiterType::Superbee, theta::Real)
    return max(0, min(1, 2*theta), min(2, theta))
end

function limiter(limiterType::Minmod, theta::Real)
    if  theta > 0
        return min(1, abs(theta))
    else
        return 0
    end
end

end
