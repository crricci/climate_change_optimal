

function optimPlanner(p; quiet = true, showObj = false)
    
    plannerProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerProblem, "tol", 1e-32)
    set_optimizer_attribute(plannerProblem, "acceptable_tol",  1e-32)
    set_optimizer_attribute(plannerProblem, "max_iter", Int(1e6))

    JuMP.@variables(plannerProblem,begin
        C1 ≥ 0.0
        C2 ≥ 0.0
        B1 ≥ 0.0
        B2 ≥ 0.0
        K1 ≥ 0.0
        K2 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(plannerProblem, begin
        gRb, (p.g∞ * Rb + p.g0) / (Rb + 1)
        hRa, (p.h∞ * Ra + p.h0) / (Ra + 1)
        fK1, K1^p.α
        fK2, K2^p.α
        D1, B1^p.θ1
        D2, gRb * B2^p.θ2
        K, K1 + K2
        D, D1 + D2
        GID, p.ηK * K - p.ηB * D
    end)
    
    # utility
    if p.σ1 == 1
        @NLexpression(plannerProblem, u1C1, log(C1))
    else
        @NLexpression(plannerProblem, u1C1, C1^(1-p.σ1)/(1-p.σ1))
    end
    if p.σ2 == 1
        @NLexpression(plannerProblem, u2C2, log(C2))
    else
        @NLexpression(plannerProblem, u2C2, C2^(1-p.σ2)/(1-p.σ2))
    end

    # budget constraint
    @NLconstraint(plannerProblem,BudgetConstraint,
    C1 + C2 + B1 + B2 + K1 + K2 + Ra + Rb ≤ p.A̅ * (fK1 + hRa * fK2))

    # objective
    @NLobjective(plannerProblem, Max, u1C1 + u2C2 - (p.γ1 + p.γ2) * p.Φ * GID)

    # verbose
    if quiet == false
        unset_silent(plannerProblem)
    else
        set_silent(plannerProblem)
    end

    JuMP.optimize!(plannerProblem)
    
    if termination_status(plannerProblem) ∉ [LOCALLY_SOLVED, OPTIMAL]
        printstyled("plannerProblem, termination status -> " * string(termination_status(plannerProblem)) * "\n",color=:red)
    end

    solValues = value.([C1,C2,B1,B2,K1,K2,Ra,Rb])
    solValues[solValues .< 0] .= 0
    solDict = Dict(zip(varNames,solValues))

    if showObj println("Objective: ",computeWelFarePlanner(p,solDict)) end

    return solDict
end

function optimPlannerExplicit(p; quiet = false)

    #check that we are in the case where there is an explicit solution 
    @assert p.α == 1

    function computeRb(p; quiet = quiet)
        RbProblem = Model(Ipopt.Optimizer)
        set_optimizer_attribute(RbProblem, "tol", 1e-16)
        JuMP.@variable(RbProblem, Rb ≥ 0)

        @NLexpression(RbProblem, gRb, (p.g∞ * Rb + p.g0) / (Rb + 1))
        @NLobjective(RbProblem, Max, -Rb + gRb^(1/(1-p.θ2)) * (p.θ2 / (p.ηK / ((p.A̅ - 1) * p.ηB)))^(1/(1-p.θ2)) * (1/p.θ2 - 1) )

        # verbose
        if quiet == false
            unset_silent(RbProblem)
        else
            set_silent(RbProblem)
        end
        JuMP.optimize!(RbProblem)
        Rb = value(Rb); Rb = Rb < 0 ? 0 : Rb
        return Rb
    end

    g(x) = (p.g∞ * x + p.g0) / (x+1)

    C1 = ( (p.A̅-1)/((p.γ1+p.γ2)*p.Φ*p.ηK) )^(1/p.σ1)
    C2 = ( (p.A̅-1)/((p.γ1+p.γ2)*p.Φ*p.ηK) )^(1/p.σ2)
    Ra = 0
    if p.ηK / p.ηB >= p.g0^p.θ2 * p.gPrime0^(1-p.θ2) * p.θ2^p.θ2 * (p.A̅ - 1)
        println("Explicit condition for Rb = 0 matched")
        Rb = 0 
    else
        Rb = computeRb(p; quiet = quiet)
    end
    B1 = (p.ηB/p.ηK * p.θ1 * (p.A̅ - 1))^(1/(1-p.θ1))
    B2 = (p.ηB/p.ηK * p.θ2 * g(Rb) * (p.A̅ - 1))^(1/(1-p.θ2))
    K1 = 1/(p.A̅ - 1) * (C1 + C2 + B1 + B2 + Rb)
    K2 = 0

    solValues = [C1,C2,B1,B2,K1,K2,Ra,Rb]
        solDict = Dict(zip(varNames,solValues))

    return solDict

end

function computeGPlanner(p)

    solDict = optimPlanner(p; quiet = true, showObj = false)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [solDict[name] for name in varNames]
    return G(K1,K2,B1,B2,Rb,p)
end

function computeWelFarePlanner(p,sol)

    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]

    # objective value
    u1C1 = p.σ1 == 1 ? log(C1) : C1^(1-p.σ1)/(1-p.σ1)
    u2C2 = p.σ2 == 1 ? log(C2) : C2^(1-p.σ2)/(1-p.σ2)
    gRb = (p.g∞ * Rb + p.g0) / (Rb + 1)
    hRa = (p.h∞ * Ra + p.h0) / (Ra + 1)
    D1 = B1^p.θ1
    D2 = gRb * B2^p.θ2
    K = K1 + K2
    D = D1 + D2
    GID = p.ηK * K - p.ηB * D

    obj = u1C1 + u2C2 - (p.γ1+p.γ2)*p.Φ*GID

    return obj
end

function computeWelFare1Planner(p,sol)

    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]

    # objective value
    u1C1 = p.σ1 == 1 ? log(C1) : C1^(1-p.σ1)/(1-p.σ1)
    gRb = (p.g∞ * Rb + p.g0) / (Rb + 1)
    D1 = B1^p.θ1
    D2 = gRb * B2^p.θ2
    K = K1 + K2
    D = D1 + D2
    GID = p.ηK * K - p.ηB * D

    obj = u1C1 - p.γ1 * p.Φ * GID

    return obj
end

function computeWelFare2Planner(p,sol)

    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]

    # objective value
    u2C2 = p.σ2 == 1 ? log(C2) : C2^(1-p.σ2)/(1-p.σ2)
    gRb = (p.g∞ * Rb + p.g0) / (Rb + 1)
    D1 = B1^p.θ1
    D2 = gRb * B2^p.θ2
    K = K1 + K2
    D = D1 + D2
    GID = p.ηK * K - p.ηB * D

    obj = u2C2 - p.γ2 * p.Φ * GID

    return obj
end



function computeAllPlanner()

    sol = optimPlanner(SOpG,quiet=true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    GComputed = computeG(SOpG,sol)

    solODE = solveODE(SOpG,DpG,GComputed)
    P,T = solODE.u[end]
    TempFinal = ComputeTemperature(DpG,P,T)
    TempInitial = ComputeTemperature(DpG,DpG.P0,DpG.T0)

    ΔP = P - DpG.P0
    ΔT = T - DpG.T0
    ΔTemp = TempFinal - TempInitial

    WelFare1 = computeWelFare1Planner(SOpG,sol)
    WelFare2 = computeWelFare2Planner(SOpG,sol)
    WelFare = computeWelFarePlanner(SOpG,sol)
    Y1 = SOpG.A̅ * K1^SOpG.α
    Y2 = SOpG.A̅ * h(Ra,SOpG) * K2^SOpG.α


    return C1,C2,B1,B2,K1,K2,Ra,Rb, GComputed, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end

function computeAllPlannerExplicit()

    sol = optimPlannerExplicit(SOpG,quiet=true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    GComputed = computeG(SOpG,sol)
    G1 = computeG1(SOpG,sol)
    G2 = computeG2(SOpG,sol)

    solODE = solveODE(SOpG,DpG,GComputed)
    P,T = solODE.u[end]
    TempFinal = ComputeTemperature(DpG,P,T)
    TempInitial = ComputeTemperature(DpG,DpG.P0,DpG.T0)

    ΔP = P - DpG.P0
    ΔT = T - DpG.T0
    ΔTemp = TempFinal - TempInitial

    WelFare1 = computeWelFare1Planner(SOpG,sol)
    WelFare2 = computeWelFare2Planner(SOpG,sol)
    WelFare = computeWelFarePlanner(SOpG,sol)
    Y1 = SOpG.A̅ * K1^SOpG.α
    Y2 = SOpG.A̅ * h(Ra,SOpG) * K2^SOpG.α
    gRb = g(Rb,SOpG)
    hRa = h(Ra,SOpG)

    return C1,C2,B1,B2,K1,K2,Ra,Rb,gRb,hRa,GComputed,G1,G2, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end