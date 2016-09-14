module asnm3

using ..Hyperbolic
using PyPlot

function exact_solution_step(x::Vector, loc::Real)
    return Vector{Float64}(x .< loc)
end


function exact_solution_cossin(x::Vector, loc::Real)
    sol = zeros(length(x))
    for i in eachindex(x)
        if x[i] >= loc && x[i] <= loc+2
            sol[i] = (1-cos(pi*(x[i]-1)))*sin(10*pi*x[i])
        end
    end
    #return sol
    return exact_solution_cos(x, loc).*sin(10*pi*x)
end

function exact_solution_cos(x::Vector, loc::Real)
    return Vector{Float64}( (x .> loc) .* (x .< loc +2)) .* (1-cos(pi*(x-1)))
end

function problem_1_a_solution()
    nCells = 200
    dt = 0.5*10/nCells
    mesh = Mesh1d(nCells, -5, 5, dt)

    println("The step-size is $(mesh.h).")
    println("The time-step is $(dt).")
    println("The CFL is $(dt/mesh.h).")

    # Initial conditions
    ICs = exact_solution_step(mesh.xCells, 0)
    print(size(mesh.xCells))
    #plot(mesh.xCells, ICs)

    # Define problems
    problemGodunov = Problem(mesh, Godunov(FluxBurgers(), Upwind()), ICs, 5, 9)
    problemLW = Problem(mesh, Godunov(FluxBurgers(), LaxWendroff()), ICs, 5, 9)
    problemMinmod = Problem(mesh, Godunov(FluxBurgers(), Minmod()), ICs, 5, 9)
    problemSuperbee = Problem(mesh, Godunov(FluxBurgers(), Superbee()), ICs, 5, 9)

    # Solve
    solve(problemGodunov)
    solve(problemLW)
    solve(problemMinmod)
    solve(problemSuperbee)

    println("The final time is $(problemGodunov.times[end])")

    # Plot
    figure(figsize=(6,3))
    subplot(1,2,1)
    plot(mesh.xCells, problemGodunov.solution[:, end], "-bo", label="Godunov")
    plot(mesh.xCells, problemLW.solution[:, end], "-yo", label="Lax-Wendroff")
    plot(mesh.xCells, problemMinmod.solution[:, end], "-ro", label="Minmod")
    plot(mesh.xCells, problemSuperbee.solution[:, end], "-go", label="Superbee")
    xlim(2.1, 2.7)
    ylim(0, 1.3)
    ylabel("u")
    xlabel("x")
    legend(loc=3, fontsize=8)
    grid()
    subplot(1,2,2)
    plot(mesh.xCells, problemGodunov.solution[:, end], "-bo", label="Godunov")
    plot(mesh.xCells, problemLW.solution[:, end], "-yo", label="Lax-Wendroff")
    plot(mesh.xCells, problemMinmod.solution[:, end], "-ro", label="Minmod")
    plot(mesh.xCells, problemSuperbee.solution[:, end], "-go", label="Superbee")
    grid()
    xlim(2.45, 2.55)
    ylim(0, 1.3)
    xlabel("x")
    tight_layout()

end


function problem_1_b_solution()
    nCells = 200
    dt = 0.5*10/nCells/2
    mesh = Mesh1d(nCells, 0, 10, dt)

    println("The step-size is $(mesh.h).")
    println("The time-step is $(dt).")
    println("The CFL is $(2*dt/mesh.h).")

    # Initial conditions
    ICs = exact_solution_cos(mesh.xCells, 1)
    #plot(mesh.xCells, ICs)

    # Define problems
    problemGodunov = Problem(mesh, Godunov(FluxLinear(1), Upwind()), ICs, 2.5, 9)
    problemLW = Problem(mesh, Godunov(FluxLinear(1), LaxWendroff()), ICs, 2.5, 9)
    problemMinmod = Problem(mesh, Godunov(FluxLinear(1), Minmod()), ICs, 2.5, 9)
    problemSuperbee = Problem(mesh, Godunov(FluxLinear(1), Superbee()), ICs, 2.5, 9)

    # Solve
    solve(problemGodunov)
    solve(problemLW)
    solve(problemMinmod)
    solve(problemSuperbee)

    println("The final time is $(problemGodunov.times[end])")

    # Plot
    figure(figsize=(6,3))
    subplot(1,2,1)
    plot(mesh.xCells, problemGodunov.solution[:, end], "-b", label="Godunov")
    plot(mesh.xCells, problemLW.solution[:, end], "-y", label="Lax-Wendroff")
    plot(mesh.xCells, problemMinmod.solution[:, end], "-r", label="Minmod")
    plot(mesh.xCells, problemSuperbee.solution[:, end], "-g", label="Superbee")
    xlim(3, 6)
    ylim(-0.1, 2.1)
    ylabel("u")
    xlabel("x")
    grid()
    subplot(1,2,2)
    plot(mesh.xCells, problemGodunov.solution[:, end], "-bo", label="Godunov")
    plot(mesh.xCells, problemLW.solution[:, end], "-yo", label="Lax-Wendroff")
    plot(mesh.xCells, problemMinmod.solution[:, end], "-ro", label="Minmod")
    plot(mesh.xCells, problemSuperbee.solution[:, end], "-go", label="Superbee")
    grid()
    xlim(4.4, 4.6)
    ylim(1.3, 2.1)
    xlabel("x")
    legend(loc=3, fontsize=8)
    tight_layout()
end

function problem_1_b_convergence()
    errorGodunovOne = Vector(0)
    errorLWOne = Vector(0)
    errorSuperbeeOne = Vector(0)
    errorMinmodOne = Vector(0)

    errorGodunovTwo = Vector(0)
    errorLWTwo = Vector(0)
    errorSuperbeeTwo = Vector(0)
    errorMinmodTwo = Vector(0)

    errorGodunovInf = Vector(0)
    errorLWInf = Vector(0)
    errorSuperbeeInf = Vector(0)
    errorMinmodInf = Vector(0)

    stepSizes = Vector(0)
    nCells = 500
    dt = 0.5*10/nCells/2
    while nCells <= 6000
        mesh = Mesh1d(nCells, 0, 10, dt)
        println("Number of cells is $(nCells), timestep $(dt).")
        println("CFL = $(2*dt/mesh.h)")
        push!(stepSizes, mesh.h)


        # Initial conditions
        ICs = exact_solution_cos(mesh.xCells, 1)

        # Define problems
        problemGodunov = Problem(mesh, Godunov(FluxLinear(1.0), Upwind()), ICs, 2.5, 9)
        problemLW = Problem(mesh, Godunov(FluxLinear(1.0), LaxWendroff()), ICs, 2.5, 9)
        problemSuperbee = Problem(mesh, Godunov(FluxLinear(1.0), Superbee()), ICs, 2.5, 9)
        problemMinmod = Problem(mesh, Godunov(FluxLinear(1.0), Minmod()), ICs, 2.5, 9)

        # Solve
        solve(problemGodunov)
        solve(problemLW)
        solve(problemSuperbee)
        solve(problemMinmod)
        exact = exact_solution_cos(mesh.xCells-problemLW.times[end], 1.)

        println("Final time is $(problemLW.times[end])")
        EGodunov = abs(exact - problemGodunov.solution[:, end])
        ELW = abs(exact - problemLW.solution[:, end])
        EMinmod = abs(exact - problemMinmod.solution[:, end])
        ESuperbee= abs(exact - problemSuperbee.solution[:, end])

        push!(errorGodunovOne, norm_p(EGodunov, 1, mesh.h))
        push!(errorLWOne, norm_p(ELW, 1, mesh.h))
        push!(errorMinmodOne, norm_p(EMinmod, 1, mesh.h))
        push!(errorSuperbeeOne, norm_p(ESuperbee, 1, mesh.h))

        push!(errorGodunovTwo, norm_p(EGodunov, 2, mesh.h))
        push!(errorLWTwo, norm_p(ELW, 2, mesh.h))
        push!(errorMinmodTwo, norm_p(EMinmod, 2, mesh.h))
        push!(errorSuperbeeTwo, norm_p(ESuperbee, 2, mesh.h))

        push!(errorGodunovInf, norm_p(EGodunov, Inf, mesh.h))
        push!(errorLWInf, norm_p(ELW, Inf, mesh.h))
        push!(errorMinmodInf, norm_p(EMinmod, Inf, mesh.h))
        push!(errorSuperbeeInf, norm_p(ESuperbee, Inf, mesh.h))

        dt /= 2
        nCells *= 2
    end

    figure(figsize=(6,6))
    subplot(1, 3, 1)
    loglog(stepSizes, 10*stepSizes, "-.k", label="Linear")
    loglog(stepSizes, 10*stepSizes.^2, "--k", label="Quadratic")
    loglog(stepSizes, errorGodunovOne, "-bo", label="Godunov")
    loglog(stepSizes, errorLWOne, "-yo", label="LW")
    loglog(stepSizes, errorMinmodOne, "-ro", label="Minmod")
    loglog(stepSizes, errorSuperbeeOne, "-go", label="Superbee")
    grid()
    legend(loc=4, fontsize=8)
    ylabel("Error")
    title("1-norm")
    xlim(0.001, 0.1)
    xlabel("h")

    subplot(1, 3, 2)
    loglog(stepSizes, 10*stepSizes, "-.k", label="Linear")
    loglog(stepSizes, 10*stepSizes.^2, "--k", label="Quadratic")
    loglog(stepSizes, errorGodunovTwo, "-bo", label="Godunov")
    loglog(stepSizes, errorLWTwo, "-yo", label="LW")
    loglog(stepSizes, errorMinmodTwo, "-ro", label="Minmod")
    loglog(stepSizes, errorSuperbeeTwo, "-go", label="Superbee")
    title("2-norm")
    grid()
    xlim(0.001, 0.1)
    xlabel("h")

    subplot(1, 3, 3)
    loglog(stepSizes, 10*stepSizes, "-.k", label="Linear")
    loglog(stepSizes, 10*stepSizes.^2, "--k", label="Quadratic")
    loglog(stepSizes, errorGodunovInf, "-bo", label="Godunov")
    loglog(stepSizes, errorLWInf, "-yo", label="LW")
    loglog(stepSizes, errorMinmodInf, "-ro", label="Minmod")
    loglog(stepSizes, errorSuperbeeInf, "-go", label="Superbee")
    title("Inf-norm")
    grid()
    xlim(0.001, 0.1)
    xlabel("h")
    tight_layout()

    dx = log(stepSizes[1]) - log(stepSizes[end])
    println("Slopes in 1-norm:")
    println("Godunov: $((log(errorGodunovOne[1]) - log(errorGodunovOne[end]))/dx)")
    println("LW: $((log(errorLWOne[1]) - log(errorLWOne[end]))/dx)")
    println("Minmod: $((log(errorMinmodOne[1]) - log(errorMinmodOne[end]))/dx)")
    println("Superbee: $((log(errorSuperbeeOne[1]) - log(errorSuperbeeOne[end]))/dx)")
    println("Slopes in 2-norm:")
    println("Godunov: $((log(errorGodunovTwo[1]) - log(errorGodunovTwo[end]))/dx)")
    println("LW: $((log(errorLWTwo[1]) - log(errorLWTwo[end]))/dx)")
    println("Minmod: $((log(errorMinmodTwo[1]) - log(errorMinmodTwo[end]))/dx)")
    println("Superbee: $((log(errorSuperbeeTwo[1]) - log(errorSuperbeeTwo[end]))/dx)")
    println("Slopes in inf-norm:")
    println("Godunov: $((log(errorGodunovInf[1]) - log(errorGodunovInf[end]))/dx)")
    println("LW: $((log(errorLWInf[1]) - log(errorLWInf[end]))/dx)")
    println("Minmod: $((log(errorMinmodInf[1]) - log(errorMinmodInf[end]))/dx)")
    println("Superbee: $((log(errorSuperbeeInf[1]) - log(errorSuperbeeInf[end]))/dx)")
end


function problem_1_a_convergence()
    errorGodunovOne = Vector(0)
    errorLWOne = Vector(0)
    errorSuperbeeOne = Vector(0)
    errorMinmodOne = Vector(0)

    errorGodunovTwo = Vector(0)
    errorLWTwo = Vector(0)
    errorSuperbeeTwo = Vector(0)
    errorMinmodTwo = Vector(0)

    errorGodunovInf = Vector(0)
    errorLWInf = Vector(0)
    errorSuperbeeInf = Vector(0)
    errorMinmodInf = Vector(0)

    stepSizes = Vector(0)
    nCells = 500
    dt = 0.5*10/nCells
    while nCells <= 4000
        mesh = Mesh1d(nCells, -5, 5, dt)
        println("Number of cells is $(nCells), timestep $(dt).")
        println("CFL = $(dt/mesh.h)")
        push!(stepSizes, mesh.h)


        # Initial conditions
        ICs = exact_solution_step(mesh.xCells, 0)

        # Define problems
        problemGodunov = Problem(mesh, Godunov(FluxBurgers(), Upwind()), ICs, 5, 9)
        problemLW = Problem(mesh, Godunov(FluxBurgers(), LaxWendroff()), ICs, 5, 9)
        problemSuperbee = Problem(mesh, Godunov(FluxBurgers(), Superbee()), ICs, 5, 9)
        problemMinmod = Problem(mesh, Godunov(FluxBurgers(), Minmod()), ICs, 5, 9)

        # Solve
        solve(problemGodunov)
        solve(problemLW)
        solve(problemSuperbee)
        solve(problemMinmod)
        exact = exact_solution_step(mesh.xCells, 2.5)

        println("Final time is $(problemLW.times[end])")
        EGodunov = abs(exact - problemGodunov.solution[:, end])
        ELW = abs(exact - problemLW.solution[:, end])
        EMinmod = abs(exact - problemMinmod.solution[:, end])
        ESuperbee= abs(exact - problemSuperbee.solution[:, end])

        push!(errorGodunovOne, norm_p(EGodunov, 1, mesh.h))
        push!(errorLWOne, norm_p(ELW, 1, mesh.h))
        push!(errorMinmodOne, norm_p(EMinmod, 1, mesh.h))
        push!(errorSuperbeeOne, norm_p(ESuperbee, 1, mesh.h))

        push!(errorGodunovTwo, norm_p(EGodunov, 2, mesh.h))
        println("Diff $((mesh.h*sum(abs(EGodunov).^2))^(1.0/2) - sqrt(mesh.h)*norm(EGodunov,2))")
        #(h*sum(abs(E).^p))^(1.0/p)
        #push!(errorGodunovTwo, sqrt(mesh.h)*norm(EGodunov,2))
        push!(errorLWTwo, norm_p(ELW, 2, mesh.h))
        push!(errorMinmodTwo, norm_p(EMinmod, 2, mesh.h))
        push!(errorSuperbeeTwo, norm_p(ESuperbee, 2, mesh.h))

        push!(errorGodunovInf, norm_p(EGodunov, Inf, mesh.h))
        push!(errorLWInf, norm_p(ELW, Inf, mesh.h))
        push!(errorMinmodInf, norm_p(EMinmod, Inf, mesh.h))
        push!(errorSuperbeeInf, norm_p(ESuperbee, Inf, mesh.h))

        #plot(mesh.xCells, problemGodunov.solution[:, end])
        #plot(mesh.xCells, EMinmod)
        #xlim(2.3, 2.6)
        dt /= 2
        nCells *= 2
    end

    figure(figsize=(6, 4.5))
    subplot(1, 3, 1)
    loglog(stepSizes, 10*stepSizes, "-.k", label="Linear")
    loglog(stepSizes, 10*stepSizes.^2, "--k", label="Quadratic")
    loglog(stepSizes, errorGodunovOne, "-bo", label="Godunov")
    loglog(stepSizes, errorLWOne, "-yo", label="LW")
    loglog(stepSizes, errorMinmodOne, "-ro", label="Minmod")
    loglog(stepSizes, errorSuperbeeOne, "-go", label="Superbee")
    grid()
    legend(loc=4, fontsize=8)
    ylabel("Error")
    title("1-norm")
    xlim(0.001, 0.1)
    xlabel("h")

    subplot(1, 3, 2)
    loglog(stepSizes, 10*stepSizes, "-.k", label="Linear")
    loglog(stepSizes, 10*stepSizes.^2, "--k", label="Quadratic")
    loglog(stepSizes, errorGodunovTwo, "-bo", label="Godunov")
    loglog(stepSizes, errorLWTwo, "-yo", label="LW")
    loglog(stepSizes, errorMinmodTwo, "-ro", label="Minmod")
    loglog(stepSizes, errorSuperbeeTwo, "-go", label="Superbee")
    title("2-norm")
    grid()
    xlim(0.001, 0.1)
    xlabel("h")

    subplot(1, 3, 3)
    loglog(stepSizes, 10*stepSizes, "-.k", label="Linear")
    loglog(stepSizes, 10*stepSizes.^2, "--k", label="Quadratic")
    loglog(stepSizes, errorGodunovInf, "-bo", label="Godunov")
    loglog(stepSizes, errorLWInf, "-yo", label="LW")
    loglog(stepSizes, errorMinmodInf, "-ro", label="Minmod")
    loglog(stepSizes, errorSuperbeeInf, "-go", label="Superbee")
    title("Inf-norm")
    grid()
    xlim(0.001, 0.1)
    xlabel("h")
    tight_layout()

    dx = log(stepSizes[1]) - log(stepSizes[end])
    println("Slopes in 1-norm:")
    println("Godunov: $((log(errorGodunovOne[1]) - log(errorGodunovOne[end]))/dx)")
    println("LW: $((log(errorLWOne[1]) - log(errorLWOne[end]))/dx)")
    println("Minmod: $((log(errorMinmodOne[1]) - log(errorMinmodOne[end]))/dx)")
    println("Superbee: $((log(errorSuperbeeOne[1]) - log(errorSuperbeeOne[end]))/dx)")
    println("Slopes in 2-norm:")
    println("Godunov: $((log(errorGodunovTwo[1]) - log(errorGodunovTwo[end]))/dx)")
    println("LW: $((log(errorLWTwo[1]) - log(errorLWTwo[end]))/dx)")
    println("Minmod: $((log(errorMinmodTwo[1]) - log(errorMinmodTwo[end]))/dx)")
    println("Superbee: $((log(errorSuperbeeTwo[1]) - log(errorSuperbeeTwo[end]))/dx)")
    println("Slopes in inf-norm:")
    println("Godunov: $((log(errorGodunovInf[1]) - log(errorGodunovInf[end]))/dx)")
    println("LW: $((log(errorLWInf[1]) - log(errorLWInf[end]))/dx)")
    println("Minmod: $((log(errorMinmodInf[1]) - log(errorMinmodInf[end]))/dx)")
    println("Superbee: $((log(errorSuperbeeInf[1]) - log(errorSuperbeeInf[end]))/dx)")
end

function problem_2()
    nCells = 1000
    dt = 0.5*10/nCells
    mesh = Mesh1d(nCells, 0, 10, dt)

    println("The step-size is $(mesh.h).")
    println("The time-step is $(dt).")
    println("The CFL is $(1*dt/mesh.h).")

    # Initial conditions
    ICs = exact_solution_cossin(mesh.xCells, 1)
    #plot(mesh.xCells, ICs)

    # Define problems
    problemGodunov = Problem(mesh, Godunov(FluxLinear(1), Upwind()), ICs, 2.5, 9)
    problemMinmod = Problem(mesh, Godunov(FluxLinear(1), Minmod()), ICs, 2.5, 9)
    problemSuperbee = Problem(mesh, Godunov(FluxLinear(1), Superbee()), ICs, 2.5, 9)

    # Solve
    solve(problemGodunov)
    solve(problemMinmod)
    solve(problemSuperbee)

    println("The final time is $(problemMinmod.times[end])")

    exact = exact_solution_cossin(mesh.xCells-problemMinmod.times[end], 1.)

    # Plot
    figure(figsize=(6,3))

    #plot(mesh.xCells, problemGodunov.solution[:, end], "-bo", label="Godunov")
    plot(mesh.xCells, problemMinmod.solution[:, end], "-r", label="Minmod")
    plot(mesh.xCells, problemSuperbee.solution[:, end], "-g", label="Superbee")
    plot(mesh.xCells, exact, "k", label="Exact")
    xlim(3.5, 5.5)
    #ylim(-0.1, 2.1)
    ylabel("u")
    xlabel("x")
    grid()
    legend(loc=3, fontsize=8)
    tight_layout()
end

#problem_1_a_convergence()
#problem_1_a_solution()
problem_2()
show()
end
