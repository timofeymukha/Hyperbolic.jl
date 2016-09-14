module Hyperbolic

using Plots
using LaTeXStrings
using Distributions

include("MeshM.jl")
include("FluxM.jl")
include("LimiterM.jl")
include("SchemeM.jl")
include("ProblemM.jl")

export Problem, solve, Scheme, LaxFriedrichs, Godunov, Roe, UpwindSimple,
numerical_flux,
flux_correction, accepts_limiter, Flux, FluxBurgers, FluxLinear,
FluxBuckleyLeverett, flux, flux_deriv, Mesh1d, Limiter, Upwind, LaxWendroff,
BeamWarming, Fromm, VanLeer, MC, Superbee, Minmod, limiter, norm_p,
exact_solution_step, exact_solution_cos, exact_solution_cossin, UpwindSimple


using .MeshM
using .FluxM
using .LimiterM
using .SchemeM
using .ProblemM

function norm_p(E::Vector{Float64}, p::Real, h::Float64)
    if p == Inf
        return norm(E, Inf)
    else
        return h^(1/p)*norm(E, p)
    end
end

gui()

end
