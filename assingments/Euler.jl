module Euler

#using PyPlot
using Plots
using DataFrames
using LaTeXStrings

using Hyperbolic.Limiters

abstract Scheme

type RoeEuler <: Scheme
    limiter::Limiter
end

type Mesh1d
    nCells::Int64
    h::Float64
    dt::Float64
    xCells::Vector{Float64}
    xFaces::Vector{Float64}
    function Mesh1d(nCells::Int64, xStart::Real, xEnd::Real, dt::Real)
        if xEnd < xStart
            error("Ill-defined interval!")
        else
            xFaces = linspace(xStart, xEnd, nCells+1)
            h = step(xFaces)
            xCells = (xFaces + 0.5*h)[1:end-1]
            return new(nCells, h, dt, xCells, xFaces)
        end
    end
end

type Problem
    mesh::Mesh1d
    method::Scheme
    ICs::Array{Float64, 2}
    solution::Array{Float64, 3}
    timeLength::Real
    times::Vector{Float64}
    timePrecision::Int

    function Problem(mesh::Mesh1d,
                     method::Scheme,
                     ICs::Array{Float64, 2},
                     timeLength::Real,
                     timePrecision::Int)
        times = Vector(1)
        times[end] = 0.0
        while times[end] < round(timeLength, timePrecision)
            push!(times, round(times[end]+mesh.dt, timePrecision))
        end


        solution = zeros(Float64, size(ICs)[1], 3, length(times))
        solution[:, :, 1] = ICs
        return new(mesh, method, ICs, solution, timeLength, times)

    end
end

function solve(problem::Problem)
    mesh = problem.mesh
    solution = problem.solution
    method = problem.method
    dt = mesh.dt
    h = mesh.h

    for i in 1:size(solution, 3)-1
        APlus, AMinus, F = numerical_flux(mesh, method, solution[:, :, i])
        solution[:, :, i+1] = solution[:, :, i] -
                          dt/h*(APlus[1:end-1, :] + AMinus[2:end, :]) -
                          dt/h*(F[2:end, :] - F[1:end-1, :])
    end

end

function numerical_flux(mesh::Mesh1d, scheme::RoeEuler, Q)
    û = zeros(size(Q, 1)+1)
    Ĥ = zeros(size(Q, 1)+1)

    APlus = zeros(size(Q, 1)+1, 3)
    AMinus = zeros(size(Q, 1)+1, 3)
    FTilde = zeros(size(Q, 1)+1, 3)
    δ = zeros(size(Q, 1)+1, 3)
    α = zeros(size(Q, 1)+1, 3)
    λ = zeros(size(Q, 1)+1, 3)
    r1 = zeros(size(Q, 1)+1, 3)
    r2 = zeros(size(Q, 1)+1, 3)
    r3 = zeros(size(Q, 1)+1, 3)

    γ = 1.4

    ρ = Q[:, 1]
    u = Q[:, 2]./ρ
    E = Q[:, 3]
    p = (γ - 1)*(E - 0.5*ρ.*u.^2)
    H = (E + p)./ρ

    û[2:end-1] = sqrt(ρ[1:end-1]).*u[1:end-1]
    û[2:end-1] += sqrt(ρ[2:end]).*u[2:end]
    û[2:end-1] ./= sqrt(ρ[1:end-1]) + sqrt(ρ[2:end])

    Ĥ[2:end-1] = sqrt(ρ[1:end-1]).*H[1:end-1] + sqrt(ρ[2:end]).*H[2:end]
    Ĥ[2:end-1] ./= sqrt(ρ[1:end-1]) + sqrt(ρ[2:end])

    ĉ = sqrt((γ - 1)*(Ĥ - 0.5*û.^2))

    δ[2:end-1, :] = Q[2:end, :] - Q[1:end-1, :]

    α[:, 2] = (γ - 1)*((Ĥ - û.^2).*δ[:, 1] + û.*δ[:, 2] - δ[:, 3])./ĉ.^2
    α[:, 3] = (δ[:, 2] + (ĉ - û).*δ[:, 1] - ĉ.*α[:, 2])./(2ĉ)
    α[:, 1] = δ[:, 1] - α[:, 2] - α[:, 3]

    α[1, :] = 0
    α[end, :] = 0

    λ[:, 1] = û - ĉ
    λ[:, 2] = û
    λ[:, 3] = û + ĉ

    if (maximum(λ)*mesh.dt/mesh.h) > 1
        println("CFL > 1 !!!")
    end

    r1[:, 1] = 1.
    r1[:, 2] = û - ĉ
    r1[:, 3] = Ĥ - û.*ĉ

    r2[:, 1] = 1.
    r2[:, 2] = û
    r2[:, 3] = 0.5*û.^2

    r3[:, 1] = 1.
    r3[:, 2] = û + ĉ
    r3[:, 3] = Ĥ + û.*ĉ

    for i in 1:size(APlus, 1)
        APlus[i, :] = max(λ[i, 1], 0)*α[i, 1]*r1[i, :]
        APlus[i, :] += max(λ[i, 2], 0)*α[i, 2]*r2[i, :]
        APlus[i, :] += max(λ[i, 3], 0)*α[i, 3]*r3[i, :]

        AMinus[i, :] = min(λ[i, 1], 0)*α[i, 1]*r1[i, :]
        AMinus[i, :] += min(λ[i, 2], 0)*α[i, 2]*r2[i, :]
        AMinus[i, :] += min(λ[i, 3], 0)*α[i, 3]*r3[i, :]
    end

    FTilde = flux_correction(mesh, scheme.limiter, α, λ, r1, r2, r3)

    return APlus, AMinus, FTilde
end

function flux_correction(mesh::Mesh1d,
                         limiterType::Limiter,
                         α, λ, r1, r2, r3)
    theta = zeros(size(α))


    for i in 2:size(theta, 1) - 1
        for j in 1:size(theta, 2)
            if λ[i, j] >= 0 && α[i, j] != 0
                theta[i, j] = α[i-1, j]/α[i, j]
            elseif α[i, j] != 0
                theta[i, j] = α[i+1, j]/α[i, j]
            end
        end
    end

    theta[1, :] = theta[2, :]

    for i in 1:size(α, 1)
        for j in 1:size(α, 2)
            α[i, j] = α[i, j] * limiter(limiterType, theta[i, j])
        end
    end

    F = zeros(size(α))

    for i in 1:size(F, 1)
        F[i, :] = 0.5*(abs(λ[i, 1])*(1-mesh.dt/mesh.h*abs(λ[i, 1]))*α[i, 1]*r1[i, :] +
                       abs(λ[i, 2])*(1-mesh.dt/mesh.h*abs(λ[i, 2]))*α[i, 2]*r2[i, :] +
                       abs(λ[i, 3])*(1-mesh.dt/mesh.h*abs(λ[i, 3]))*α[i, 3]*r3[i, :])
    end

    return F
end


# Euler
# Define the time-step
dt = 0.001

# Define the mesh
mesh = Mesh1d(400, 0, 1, dt)
gamma = 1.4
x = mesh.xCells

xMid = 0.5*(x[1] + x[end])
# Define the initial conditions
ICs = zeros(mesh.nCells, 3)
ICs[:, 1] = (x .<= xMid)*1 + (x .> xMid)*0.125
ICs[:, 3] = (x .<= xMid)*1/(gamma-1) + (x .> xMid)*0.1/(gamma - 1)

problemGodunov = Problem(mesh, RoeEuler(Upwind()), ICs, 0.2, 10)
problemSuperbee = Problem(mesh, RoeEuler(Superbee()), ICs, 0.2, 10)
solve(problemGodunov)
solve(problemSuperbee)

analytical = readtable("exact.txt", separator=' ')
pyclaw = readtable("pyclaw.txt", separator= ' ')
pyclaw[2] = pyclaw[2]./pyclaw[1]

rhoGodunov = problemGodunov.solution[:, 1, :]
uGodunov = problemGodunov.solution[:, 2, :]./rhoGodunov
EGodunov = problemGodunov.solution[:, 3, :]
pGodunov = (gamma - 1)*EGodunov - 0.5*(gamma - 1).*rhoGodunov.*uGodunov.^2
eGodunov = (EGodunov - 0.5*rhoGodunov.*uGodunov.^2)./rhoGodunov

rhoSuperbee = problemSuperbee.solution[:, 1, :]
uSuperbee = problemSuperbee.solution[:, 2, :]./rhoSuperbee
ESuperbee = problemSuperbee.solution[:, 3, :]
pSuperbee = (gamma - 1)*ESuperbee - 0.5*(gamma - 1).*rhoSuperbee.*uSuperbee.^2
eSuperbee = (ESuperbee - 0.5*rhoSuperbee.*uSuperbee.^2)./rhoSuperbee

#figure(figsize=(5,4))
#subplot(2, 2, 1)
#plot(analytical_u[1], analytical_u[2])
#plot(analytical[1], analytical[4])
pyplot(size=(600,400))


#plot(Mesh1d(400, 0, 1, dt).xCells, pyclaw[2], "r")
function plot_godunov_superbee()
    plotU = plot(mesh.xCells, [ uGodunov[:, end]  uSuperbee[:, end]],
                 color=[:blue :red],label=[L"\mathrm{Godunov}" L"\mathrm{Superbee}"],
                 title=L"\mathrm{Velocity}", titlefont=font(8),
                 ylabel=L"u")

    plotP = plot(mesh.xCells, [ pGodunov[:, end]  pSuperbee[:, end]],
                 color=[:blue :red],
                 title=L"\mathrm{Pressure}", titlefont=font(8),
                 ylabel=L"p", legend=:none)

    plotRho = plot(mesh.xCells, [ rhoGodunov[:, end]  rhoSuperbee[:, end]],
                   color=[:blue :red],
                   title=L"\mathrm{Density}", titlefont=font(8),
                   xlabel=L"{x}", ylabel=L"\rho", legend=:none)
    plot(plotU, plotP, plotRho, layout=(3,1), show=true)
    gui()
end

function plot_godunov_bench()
    plotU = plot(mesh.xCells, uGodunov[:, end],
                 color=:blue, label=L"\mathrm{Godunov}",
                 title=L"\mathrm{Velocity}", titlefont=font(8),
                 ylabel=L"u")
    plot!(analytical[1], analytical[4], color=:black, label=L"\mathrm{Exact}")

    plotP = plot(mesh.xCells, pGodunov[:, end],
                 color=:blue,
                 title=L"\mathrm{Pressure}", titlefont=font(8),
                 ylabel=L"p", legend=:none)

    plot!(analytical[1], analytical[2], color=:black, label=L"\mathrm{Exact}")

    plotRho = plot(mesh.xCells, rhoGodunov[:, end],
                  color=:blue,
                  title=L"\mathrm{Density}", titlefont=font(8),
                  xlabel=L"{x}", ylabel=L"\rho", legend=:none)

    plot!(analytical[1], analytical[3], color=:black, label=L"\mathrm{Exact}")
    plot(plotU, plotP, plotRho, layout=(3,1), show=true)
    gui()
end

#plot_godunov_superbee()
plot(mesh.xCells, sqrt(1.4*pSuperbee[:,end]./rhoSuperbee[:, end]), show=true,
      label=L"c", color=:black, size=(400,200))
plot!(mesh.xCells, uSuperbee[:, end], show=true,
            label=L"u", color=:red)
xlabel!(L"x")
end
