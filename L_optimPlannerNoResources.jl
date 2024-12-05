
function optimPlannerNoResources(p; RaImpact = 1.0, RbImpact = 1.0, quiet = true, start = nothing, showObj = false)
    
    plannerNRProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerNRProblem, "tol", GLOBAL_TOL)
    set_optimizer_attribute(plannerNRProblem, "acceptable_tol",  GLOBAL_TOL)
    set_optimizer_attribute(plannerNRProblem, "max_iter", GLOBAL_MAX_IT)

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
    gRb, (p.g∞ * RbImpact*Rb + p.g0) / (RbImpact*Rb + 1)
    hRa, (p.h∞ * RaImpact*Ra + p.h0) / (RaImpact*Ra + 1)
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
    solValues[solValues .< 1e-6] .= 0
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

function optimPlannerNoResourcesRobust(SOp, Dp; quiet = true)

    SOpInner = deepcopy(SOp)

    #check that we are in the case where there is a soluton, i.e. AK log case
    @assert SOpInner.α == 1
    @assert SOpInner.σ1 == 1
    @assert SOpInner.σ2 == 1

    candidateSol = optimPlannerNoResourcesExplicit(SOpInner,quiet=true)
    γ1 = 0.0
    γ2 = 0.0
    
    err = 1.0; Nit = 0
    while (err > GLOBAL_TOL) & (Nit < 50)
        
        candidateSolPre = candidateSol
        γ1Pre = SOpInner.γ1 
        γ2Pre = SOpInner.γ2

        γ1,γ2 = bestResponseNatureNoResources(SOpInner,Dp,candidateSolPre; quiet = true)
        SOpInner.γ1 = γ1
        SOpInner.γ2 = γ2
        candidateSol = optimPlannerNoResourcesExplicit(SOpInner,quiet=true)

        err = norm([γ1,γ2,[candidateSol[name] for name in varNames]...] - [γ1Pre,γ2Pre,[candidateSolPre[name] for name in varNames]...],Inf)
        Nit = Nit + 1 

        if quiet == false @show err, Nit end
    end


    if Nit == GLOBAL_MAX_IT
        printstyled("nashProblem -> maximum iteration reached" * "\n",color=:red)
    end

    solValues = [candidateSol[name] for name in varNames]
    solDict = Dict(zip(varNames,solValues))

    return solDict, γ1, γ2

end

function bestResponseNatureNoResources(SOp,Dp,sol; quiet = false)

    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    g(x) = (SOp.g∞ * x + SOp.g0) / (x+1)
    h(x) = (SOp.h∞ * x + SOp.h0) / (x+1)

    RobustProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(RobustProblem, "tol", GLOBAL_TOL)

    JuMP.@variables(RobustProblem, begin
    γ1o ≥ 0
    γ2o ≥ 0
    end)

    @NLexpression(RobustProblem, hRa, (SOp.h∞ * Ra + SOp.h0) / (Ra+1))

    Λ1 = (Dp.S̅/SOp.ρ - Dp.P0/SOp.ρ - Dp.T0/(SOp.ρ+SOp.ϕ)) +
    - SOp.Φ*SOp.ηK/SOp.ρ * (Ra+Rb+B1+B2) + 
    + SOp.Φ*SOp.ηB/SOp.ρ * (B1^SOp.θ1 + g(Rb)*B2^SOp.θ2)

    @NLobjective(RobustProblem, Min,
    1/SOp.ρ*log( (SOp.A̅-1) / ((γ1o+γ2o)*SOp.Φ*SOp.ηK) * (hRa*SOp.A̅-1) / ((γ1o+γ2o)*SOp.Φ*SOp.ηK) ) + 
    +(γ1o + γ2o) * Λ1 + hRa*SOp.A̅ + SOp.αR/SOp.ρ*((γ1o-SOp.γ̂1)^2 + (γ2o-SOp.γ̂2)^2) )
    
    # verbose
    if quiet == false
        unset_silent(RobustProblem)
    else
        set_silent(RobustProblem)
    end
    JuMP.optimize!(RobustProblem)

    γ1, γ2 = value.([γ1o,γ2o])
    
    return γ1, γ2

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




function computeAllPlannerNoResources(SOp,Dp; RaImpact = 1.0, RbImpact = 1.0)

    sol = optimPlannerNoResources(SOp,quiet=true,RaImpact=RaImpact,RbImpact=RbImpact)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    GComputed = computeG(SOp,sol)
    G1 = computeG1(SOp,sol)
    G2 = computeG2(SOp,sol)

    solODE = solveODE(SOp,Dp,GComputed)
    P,T = solODE.u[end]
    TempFinal = ComputeTemperature(Dp,P,T)
    TempInitial = ComputeTemperature(Dp,Dp.P0,Dp.T0)

    ΔP = P - Dp.P0
    ΔT = T - Dp.T0
    ΔTemp = TempFinal - TempInitial

    WelFare1 = computeWelFare1PlannerNoResources(SOp,sol)
    WelFare2 = computeWelFare2PlannerNoResources(SOp,sol)
    WelFare = computeWelFarePlannerNoResources(SOp,sol)
    Y1 = SOp.A̅ * K1^SOp.α
    Y2 = SOp.A̅ * h(Ra,SOp) * K2^SOp.α
    gRb = g(Rb,SOp)
    hRa = h(Ra,SOp)

    return C1,C2,B1,B2,K1,K2,Ra,Rb,gRb,hRa, GComputed,G1,G2, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end

function computeAllPlannerNoResourcesExplicit(SOp,Dp)

    sol = optimPlannerNoResourcesExplicit(SOp,quiet=true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    GComputed = computeG(SOp,sol)
    G1 = computeG1(SOp,sol)
    G2 = computeG2(SOp,sol)

    solODE = solveODE(SOp,Dp,GComputed)
    P,T = solODE.u[end]
    TempFinal = ComputeTemperature(Dp,P,T)
    TempInitial = ComputeTemperature(Dp,Dp.P0,Dp.T0)

    ΔP = P - Dp.P0
    ΔT = T - Dp.T0
    ΔTemp = TempFinal - TempInitial

    WelFare1 = computeWelFare1PlannerNoResources(SOp,sol)
    WelFare2 = computeWelFare2PlannerNoResources(SOp,sol)
    WelFare = computeWelFarePlannerNoResources(SOp,sol)
    Y1 = SOp.A̅ * K1^SOp.α
    Y2 = SOp.A̅ * h(Ra,SOp) * K2^SOp.α
    gRb = g(Rb,SOp)
    hRa = h(Ra,SOp)

    return C1,C2,B1,B2,K1,K2,Ra,Rb,gRb,hRa, GComputed,G1,G2, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end
