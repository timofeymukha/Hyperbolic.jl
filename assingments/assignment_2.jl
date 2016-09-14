module asnm2

#include("..\\Hyperbolic.jl")

using ..Hyperbolic
using LaTeXStrings
using Distributions
using Plots


function problem_1_a()
    dt = 0.02
    mesh = Mesh1d(250, -5, 5, dt)
    ICs = 0.2*ones(mesh.nCells)
    ICs[1:Int(ceil(length(ICs)/2.0))] = 1
    CFL = mesh.dt/mesh.h

    problemGodunov = Problem(mesh, Godunov(FluxBurgers(), Upwind()), ICs, 5.0, 5)
    solve(problemGodunov)
    problemLW= Problem(mesh, Godunov(FluxBurgers(), LaxWendroff()), ICs, 5.0, 5)
    solve(problemLW)
    problemLF = Problem(mesh, LaxFriedrichs(FluxBurgers()), ICs, 5.0, 5)
    solve(problemLF)
    problemUpwind = Problem(mesh, UpwindSimple(FluxBurgers()), ICs, 5.0, 5)
    solve(problemUpwind)
    problemRoe = Problem(mesh, Roe(FluxBurgers()), ICs, 5.0, 5)
    solve(problemRoe)

    plot(mesh.xCells, problemGodunov.solution[:, end], m=:circle, markersize=2,
    c=:blue, label=L"\mathrm{Godunov}", size=(600,200))
    plot!(mesh.xCells, problemLF.solution[:, end], m=:circle, markersize=2,
    c=:red, label=L"$\mathrm{Lax}$-$\mathrm{Friedrichs}$")
    plot!(mesh.xCells, problemLW.solution[:, end], m=:circle, markersize=2,
    c=:green, label=L"$\mathrm{Lax}$-$\mathrm{Wendroff}$")
    plot!(mesh.xCells, problemUpwind.solution[:, end], m=:circle, markersize=2,
    c=:magenta, label=L"$\mathrm{Non}$-$\mathrm{conservative}$")
    plot!(mesh.xCells, problemRoe.solution[:, end], m=:circle, markersize=2,
    c=:yellow, label=L"\mathrm{Roe}")
    xlabel!(L"x")
    ylabel!(L"u")
    xlims!((2, 3.5))
end

function problem_1_b()
    dt = 0.02
    mesh = Mesh1d(250, -5, 5, dt)
    ICs = (mesh.xCells .> 0) * 1. + (mesh.xCells .< 0) * (-1)
    CFL = mesh.dt/mesh.h
    println(CFL)

    problemGodunov = Problem(mesh, Godunov(FluxBurgers(), Upwind()), ICs, 1.0, 5)
    solve(problemGodunov)
    problemLW= Problem(mesh, Godunov(FluxBurgers(), LaxWendroff()), ICs, 1.0, 5)
    solve(problemLW)
    problemLF = Problem(mesh, LaxFriedrichs(FluxBurgers()), ICs, 1.0, 5)
    solve(problemLF)
    problemUpwind = Problem(mesh, UpwindSimple(FluxBurgers()), ICs, 1.0, 5)
    solve(problemUpwind)
    problemRoe = Problem(mesh, Roe(FluxBurgers()), ICs, 1.0, 5)
    solve(problemRoe)

    plot(mesh.xCells, problemGodunov.solution[:, end],
    c=:blue, label=L"\mathrm{Godunov}", size=(600,200))
    plot!(mesh.xCells, problemLF.solution[:, end],
    c=:red, label=L"$\mathrm{Lax}$-$\mathrm{Friedrichs}$")
    plot!(mesh.xCells, problemLW.solution[:, end],
    c=:green, label=L"$\mathrm{Lax}$-$\mathrm{Wendroff}$")
    plot!(mesh.xCells, problemUpwind.solution[:, end],
    c=:magenta, label=L"$\mathrm{Non}$-$\mathrm{conservative}$")
    plot!(mesh.xCells, problemRoe.solution[:, end],
    c=:yellow, label=L"\mathrm{Roe}")
    xlabel!(L"x")
    ylabel!(L"u")
    xlims!((-2.5, 2.5))
end

function problem_2()
    dt = 0.001
    mesh = Mesh1d(1500, -1, 5, dt)
    ICs = (mesh.xCells .> 0) * 0. + (mesh.xCells .< 0) * (1)
    CFL = maximum(flux_deriv(FluxBuckleyLeverett(),mesh.xCells))*mesh.dt/mesh.h
    println(CFL)

    problemLW= Problem(mesh, Godunov(FluxBuckleyLeverett(), LaxWendroff()), ICs, 1.0, 10)
    solve(problemLW)
    problemLF = Problem(mesh, LaxFriedrichs(FluxBuckleyLeverett()), ICs, 1.0, 10)
    solve(problemLF)
    problemRoe = Problem(mesh, Roe(FluxBuckleyLeverett()), ICs, 1, 10)
    solve(problemRoe)

    plot(mesh.xCells, problemLF.solution[:, end], size=(600,200),
    c=:red, label=L"$\mathrm{Lax}$-$\mathrm{Friedrich}$")
    plot!(mesh.xCells, problemLW.solution[:, end],
    c=:green, label=L"$\mathrm{Lax}$-$\mathrm{Wendroff}$")
    plot!(mesh.xCells, problemRoe.solution[:, end],
    c=:yellow, label=L"\mathrm{Roe}")
    xlabel!(L"x")
    ylabel!(L"u")
    xlims!((-0.5, 2))
    #show()

    plot(mesh.xCells, (flux(FluxBuckleyLeverett(), problemRoe.solution[:, end])-1)./(problemRoe.solution[:, end] - 1))
    plot!(mesh.xCells, flux(FluxBuckleyLeverett(), problemRoe.solution[:, end])./problemRoe.solution[:, end])
    plot(mesh.xCells, (flux(FluxBuckleyLeverett(), problemLF.solution[:, end])-1)./(problemLF.solution[:, end] - 1), c=:red)
    plot!(mesh.xCells, flux(FluxBuckleyLeverett(), problemLF.solution[:, end])./problemLF.solution[:, end], c=:red, linestyle=:dash)
    xlims!((-0.5, 2))
end

function exact_solution_cos(x::Vector, loc::Real)
    return Vector{Float64}( (x .> loc) .* (x .< loc +2)) .* (1-cos(pi*(x-1)))
end

function problem_3()
    dt = 0.05
    mesh = Mesh1d(500, -50, 150, dt)
    ICs = 0.75*pdf(Normal(-25, 5), mesh.xCells)/maximum(pdf(Normal(-25, 5), mesh.xCells)) +0.25
    a = zeros(length(mesh.xCells))
    #ICs = 0.2*ones(size(mesh.xCells))

    for i in eachindex(a)
        if mesh.xCells[i] <= 0
            a[i] = 1
        elseif mesh.xCells[i] <= 50
            a[i] = 2
        elseif mesh.xCells[i] <= 100
            a[i] = 1
        #elseif mesh.xCells[i] <= 150
        #    a[i] = 2
        #elseif mesh.xCells[i] <= 200
        #    a[i] = 1
        else
            a[i] = 2
        end

    end

    problemLF = Problem(mesh, LocalLaxFriedrichs(FluxTrafficMultilane(1, a)), ICs, 100.0, 10)
    solve(problemLF)

    solPlot = plot(mesh.xCells, problemLF.solution[:, 1], linestyle=:dash,
    c=:red, label=L"$\mathrm{Local Lax}$-$\mathrm{Friedrichs}$")
    plot!(mesh.xCells, problemLF.solution[:, end],
    c=:red, label=L"$\mathrm{Local Lax}$-$\mathrm{Friedrichs}$")

    ylabel!(L"\rho")
    #xlims!((-50, 30))
    ylims!((0, 1))

    aPlot = plot(mesh.xCells, a, label=L"a")
    ylims!((0.9, 2.1))
    xlabel!(L"x")
    ylabel!(L"a")

    plot(solPlot, aPlot, layout=grid(2,1, heights=[0.8, 0.2]))
end

problem_2()

#plot(linspace(0, 1), flux(FluxBuckleyLeverett(), linspace(0,1))./linspace(0,1), legend=:none, size=(400,200))
#xlabel!(L"u^*")
#ylabel!(L"f(u^*)/u^*")

end
