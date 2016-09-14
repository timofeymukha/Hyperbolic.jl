module ProblemM

using Hyperbolic.SchemeM
using Hyperbolic.LimiterM
using Hyperbolic.MeshM

export Problem, solve

type Problem
    mesh::Mesh1d
    method::Scheme
    ICs::Array{Float64, 1}
    solution::Array{Float64, 2}
    timeLength::Real
    times::Vector{Float64}
    timePrecision::Int

    function Problem(mesh::Mesh1d,
                     method::Scheme,
                     ICs::Array{Float64, 1},
                     timeLength::Real,
                     timePrecision::Int)
        times = Vector(1)
        times[end] = 0.0
        while times[end] < round(timeLength, timePrecision)
            push!(times, round(times[end]+mesh.dt, timePrecision))
        end

        solution = zeros(Float64, length(ICs), length(times))
        solution[:, 1] = ICs

        return new(mesh, method, ICs, solution, timeLength, times)

    end
end

function solve(problem::Problem)
    mesh = problem.mesh
    solution = problem.solution
    method = problem.method
    dt = mesh.dt
    h = mesh.h

    for i in 1:size(solution, 2)-1
        APlus, AMinus, F = numerical_flux(mesh, method, solution[:, i])
        solution[:, i+1] = solution[:, i] -
                          dt/h*(APlus[1:end-1] + AMinus[2:end]) -
                          dt/h*(F[2:end] - F[1:end-1])
    end

end
end
