
function optimPlannerNoResources(p; quiet = true, start = nothing, showObj = false)
    
    plannerNRProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerNRProblem, "tol", 1e-32)
    set_optimizer_attribute(plannerNRProblem, "acceptable_tol",  1e-32)
    set_optimizer_attribute(plannerNRProblem, "max_iter", Int(1e6))

    JuMP.@variables(plannerNRProblem,begin
        C1 ≥ 0.0
        C2 ≥ 0.0
        B1 ≥ 0.0
        B2 ≥ 0.0
        K1 ≥ 0.0
        K2 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    if !isnothing(start)
        set_start_value(C1,start["C1"])
        set_start_value(C2,start["C2"])
        set_start_value(B1,start["B1"])
        set_start_value(B2,start["B2"])
        set_start_value(K1,start["K1"])
        set_start_value(K2,start["K2"])
        set_start_value(Ra,start["Ra"])
        set_start_value(Rb,start["Rb"])
    end

    # subexpressions
    @NLexpressions(plannerNRProblem, begin
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
        @NLexpression(plannerNRProblem, u1C1, log(C1))
    else
        @NLexpression(plannerNRProblem, u1C1, C1^(1-p.σ1)/(1-p.σ1))
    end
    if p.σ2 == 1
        @NLexpression(plannerNRProblem, u2C2, log(C2))
    else
        @NLexpression(plannerNRProblem, u2C2, C2^(1-p.σ2)/(1-p.σ2))
    end

    # budget constraint 1
    @NLconstraint(plannerNRProblem,BudgetConstraint1,
    C1 + B1 + K1 + Ra + Rb ≤ p.A̅ * fK1)

    # budget constraint 2
    @NLconstraint(plannerNRProblem,BudgetConstraint2,
    C2 + B2 + K2 ≤ p.A̅ * hRa * fK2)
    
    # objective
    @NLobjective(plannerNRProblem, Max, u1C1 + u2C2 - (p.γ1 + p.γ2) * p.Φ * GID)

    # verbose
    if quiet == false
        unset_silent(plannerNRProblem)
    else
        set_silent(plannerNRProblem)
    end

    JuMP.optimize!(plannerNRProblem)
        
    if termination_status(plannerNRProblem) ∉ [LOCALLY_SOLVED, OPTIMAL]
        printstyled("plannerNRProblem, termination status -> " * string(termination_status(plannerNRProblem)) * "\n",color=:red)
    end

    solValues = value.([C1,C2,B1,B2,K1,K2,Ra,Rb])
    solValues[solValues .< 0] .= 0
    solDict = Dict(zip(varNames,solValues))

    if showObj println("Objective: ",computeWelFarePlannerNoResources(p,solDict)) end
    
    return solDict
end

function optimPlannerNoResourcesExplicit(p; quiet = true)

    #check that we are in the case where there is an explicit solution 
    @assert p.α == 1

    function computeRaRb(p; quiet = quiet)
        RaRbProblem = Model(Ipopt.Optimizer)
        # set_optimizer_attribute(RaRbProblem, "tol", 1e-16)
        JuMP.@variable(RaRbProblem, Ra ≥ 0)
        JuMP.@variable(RaRbProblem, Rb ≥ 0)


        @NLexpression(RaRbProblem, gRb, (p.g∞ * Rb + p.g0) / (Rb + 1))
        @NLexpression(RaRbProblem, hRa, (p.h∞ * Ra + p.h0) / (Ra + 1))
        if p.σ2 == 1
            @NLobjective(RaRbProblem, Max, 
            1/((p.γ1+p.γ2)*p.Φ)*log((p.A̅*hRa - 1)/((p.γ1+p.γ2)*p.Φ*p.ηK)) + 
            (1-p.θ2)*p.ηB*gRb * (p.ηK/ ((p.ηB * p.θ2 * gRb) * (p.A̅*hRa - 1)))^(p.θ2/(p.θ2 - 1)) - p.ηK/(p.A̅-1)*(Ra+Rb) )
        else
            @NLobjective(RaRbProblem, Max, 
            p.σ2/(1-p.σ2) * (1/((p.γ1+p.γ2)*p.Φ*p.ηK))^(1/p.σ2) * (p.A̅*hRa - 1)^((1-p.σ2)/p.σ2) + 
            (1-p.θ2)*p.ηB*gRb * (p.ηK/ ((p.ηB * p.θ2 * gRb) * (p.A̅*hRa - 1)))^(p.θ2/(p.θ2 - 1)) -
            p.ηK/(p.A̅-1)*(Ra+Rb) )
        end

        # verbose
        if quiet == false
            unset_silent(RaRbProblem)
        else
            set_silent(RaRbProblem)
        end
        JuMP.optimize!(RaRbProblem)
        Rb = value(Rb); Rb = Rb < 0 ? 0 : Rb
        Ra = value(Ra); Ra = Ra < 0 ? 0 : Ra

        return Ra,Rb
    end

    g(x) = (p.g∞ * x + p.g0) / (x+1)
    h(x) = (p.h∞ * x + p.h0) / (x+1)
    
    Ra,Rb = computeRaRb(p; quiet=quiet)
    C1 = ( (p.γ1+p.γ2) * p.Φ * p.ηK / (p.A̅ - 1) )^(-1/p.σ1)
    C2 = ((p.A̅ * h(Ra) - 1)/((p.γ1+p.γ2) * p.Φ * p.ηK))^(1/p.σ2)
    B1 = (p.ηB/p.ηK * p.θ1 * (p.A̅ - 1))^(1/(1-p.θ1))
    B2 = (p.ηK/ ((p.ηB * p.θ2 * g(Rb)) * (p.A̅ * h(Ra) - 1)))^(1/(p.θ2 - 1))
    K1 = 1/(p.A̅ - 1) * (C1 + B1 + Ra + Rb)
    K2 = 1/(p.A̅*h(Ra) - 1) * (C2 + B2)

    
    solValues = [C1,C2,B1,B2,K1,K2,Ra,Rb]
    solDict = Dict(zip(varNames,solValues))

    # println("Objective: ", computeWelFarePlannerNoResources(p,solDict))
    
    return solDict
 
end

function computeGPlannerNoResources(p)

    solDict = optimPlannerNoResources(p; quiet = true, showObj = false)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [solDict[name] for name in varNames]
    return G(K1,K2,B1,B2,Rb,p)
end

function computeWelFarePlannerNoResources(p,sol)

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

function computeWelFare1PlannerNoResources(p,sol)

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

function computeWelFare2PlannerNoResources(p,sol)

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




function computeAllPlannerNoResources()

    sol = optimPlannerNoResources(SOpG,quiet=true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    GComputed = computeG(SOpG,sol)


    solODE = solveODE(SOpG,DpG,GComputed)
    P,T = solODE.u[end]
    TempFinal = ComputeTemperature(DpG,P,T)
    TempInitial = ComputeTemperature(DpG,DpG.P0,DpG.T0)

    ΔP = P - DpG.P0
    ΔT = T - DpG.T0
    ΔTemp = TempFinal - TempInitial

    WelFare1 = computeWelFare1PlannerNoResources(SOpG,sol)
    WelFare2 = computeWelFare2PlannerNoResources(SOpG,sol)
    WelFare = computeWelFarePlannerNoResources(SOpG,sol)
    Y1 = SOpG.A̅ * K1^SOpG.α
    Y2 = SOpG.A̅ * h(Ra,SOpG) * K2^SOpG.α


    return C1,C2,B1,B2,K1,K2,Ra,Rb, GComputed, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end

function computeAllPlannerNoResourcesExplicit()

    sol = optimPlannerNoResourcesExplicit(SOpG,quiet=true)
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

    WelFare1 = computeWelFare1PlannerNoResources(SOpG,sol)
    WelFare2 = computeWelFare2PlannerNoResources(SOpG,sol)
    WelFare = computeWelFarePlannerNoResources(SOpG,sol)
    Y1 = SOpG.A̅ * K1^SOpG.α
    Y2 = SOpG.A̅ * h(Ra,SOpG) * K2^SOpG.α
    gRb = g(Rb,SOpG)
    hRa = h(Ra,SOpG)

    return C1,C2,B1,B2,K1,K2,Ra,Rb,gRb,hRa, GComputed,G1,G2, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end
