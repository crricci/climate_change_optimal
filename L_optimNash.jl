function optimNash(p; RaImpact = 1.0, RbImpact = 1.0, TOL = GLOBAL_TOL, MAX_IT = GLOBAL_MAX_IT, quiet = true)
    # by iterated best response

    # initial point by social planner optimum
    solPlanner = optimPlanner(p; quiet = true, showObj = false)
    B2,K2 = solPlanner["B2"],solPlanner["K2"]
    
    C1,B1,K1,Ra,Rb = NashNation1(p,B2,K2,RaImpact,RbImpact)
    C2, B2, K2 = NashNation2(p,B1,K1,Ra,Rb,RaImpact,RbImpact)
    candidateSol = [C1,C2,B1,B2,K1,K2,Ra,Rb]

    err = 1.0; Nit = 0
    while (err > GLOBAL_TOL) & (Nit < 50)
        
        C1,B1,K1,Ra,Rb = NashNation1(p,B2,K2,RaImpact,RbImpact)
        C2, B2, K2 = NashNation2(p,B1,K1,Ra,Rb,RaImpact,RbImpact)
        
        err = norm(candidateSol - [C1,C2,B1,B2,K1,K2,Ra,Rb],Inf)
        Nit = Nit + 1 

        candidateSol = [C1,C2,B1,B2,K1,K2,Ra,Rb]
        
        if quiet == false 
            print("Nash Case: ")
            @show err, Nit
         end
    end
    
    if Nit == GLOBAL_MAX_IT
        printstyled("nashProblem -> maximum iteration reached" * "\n",color=:red)
    end

    solDict = Dict(zip(varNames,candidateSol))

    return solDict
end

function NashNation1(p,B2,K2,RaImpact,RbImpact)

    nash1Problem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(nash1Problem, "tol", GLOBAL_TOL)
    set_optimizer_attribute(nash1Problem, "acceptable_tol",  GLOBAL_TOL)
    set_optimizer_attribute(nash1Problem, "max_iter", GLOBAL_MAX_IT)

    JuMP.@variables(nash1Problem,begin
        C1 ≥ 0.0
        B1 ≥ 0.0
        K1 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(nash1Problem, begin
    gRb, (p.g∞ * RbImpact*Rb + p.g0) / (RbImpact*Rb + 1)
    hRa, (p.h∞ * RaImpact*Ra + p.h0) / (RaImpact*Ra + 1)
    fK1, K1^p.α
    D1, B1^p.θ1
    D2, gRb * B2^p.θ2
    K, K1 + K2
    D, D1 + D2
    GID, p.ηK * K - p.ηB * D
    end)
    
    # utility
    if p.σ1 == 1
        @NLexpression(nash1Problem, u1C1, log(C1))
    else
        @NLexpression(nash1Problem, u1C1, C1^(1-p.σ1)/(1-p.σ1))
    end

    # budget constraint
    @NLconstraint(nash1Problem,BudgetConstraint,
    C1 + B1 + K1 + Ra + Rb ≤ p.A̅ * fK1)

    # objective
    @NLobjective(nash1Problem, Max, u1C1 - p.γ1 * p.Φ * GID)

    set_silent(nash1Problem)
    JuMP.optimize!(nash1Problem)
    
    solValues = value.([C1,B1,K1,Ra,Rb])
    solValues[solValues .< 1e-6] .= 0
    C1,B1,K1,Ra,Rb = solValues

    return C1,B1,K1,Ra,Rb
end

function NashNation2(p,B1,K1,Ra,Rb,RaImpact,RbImpact)

    nash2Problem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(nash2Problem, "tol", GLOBAL_TOL)
    set_optimizer_attribute(nash2Problem, "acceptable_tol",  GLOBAL_TOL)
    set_optimizer_attribute(nash2Problem, "max_iter", GLOBAL_MAX_IT)

    JuMP.@variables(nash2Problem,begin
        C2 ≥ 0.0
        B2 ≥ 0.0
        K2 ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(nash2Problem, begin
    gRb, (p.g∞ * RbImpact*Rb + p.g0) / (RbImpact*Rb + 1)
    hRa, (p.h∞ * RaImpact*Ra + p.h0) / (RaImpact*Ra + 1)
    fK2, K2^p.α
    D1, B1^p.θ1
    D2, gRb * B2^p.θ2
    K, K1 + K2
    D, D1 + D2
    GID, p.ηK * K - p.ηB * D
    end)
    
    # utility
    if p.σ2 == 1
        @NLexpression(nash2Problem, u2C2, log(C2))
    else
        @NLexpression(nash2Problem, u2C2, C2^(1-p.σ2)/(1-p.σ2))
    end

    # budget constraint
    @NLconstraint(nash2Problem,BudgetConstraint,
    C2 + B2 + K2  ≤ p.A̅ * hRa * fK2)

    # objective
    @NLobjective(nash2Problem, Max, u2C2 -  p.γ2 * p.Φ * GID)

    set_silent(nash2Problem)
    JuMP.optimize!(nash2Problem)

    solValues = value.([C2, B2, K2])
    solValues[solValues .< 1e-6] .= 0
    C2, B2, K2 = solValues

    return C2, B2, K2
end

function optimNashRobust(SOp,Dp; quiet = true)
    SOpInner = deepcopy(SOp)

    #check that we are in the case where there is a soluton, i.e. AK log case
    @assert SOp.α == 1
    @assert SOp.σ1 == 1
    @assert SOp.σ2 == 1

    # initialization
    γ1 = SOpInner.γ1
    γ2 = SOpInner.γ2
    solNashStart = optimNashExplicit(SOpInner; quiet = true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [solNashStart[name] for name in varNames]

    err = 1.0; Nit = 0
    while (err > GLOBAL_TOL) & (Nit < 50)
        
        γ1Pre,γ2Pre = γ1,γ2 

        γ1 = NashRobustNation1gamma1(SOpInner,Dp,B1,B2,Rb,γ2)
        γ2 = NashRobusNation2gamma2(SOpInner,Dp,B1,B2,Rb,γ1)
        
        err = norm([γ1,γ2] - [γ1Pre,γ2Pre],Inf)
        Nit = Nit + 1 

        if quiet == false 
            print("Nash Case: ")
            @show err, Nit
         end
    end
    
    if Nit == GLOBAL_MAX_IT
        printstyled("nashProblem -> maximum iteration reached" * "\n",color=:red)
    end

    # after finding the correct gamma recompute the strategies of the nations 
    SOpInner.γ1 = γ1
    SOpInner.γ2 = γ2
    solNash = optimNashExplicit(SOpInner; quiet = true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [solNash[name] for name in varNames]

    solDict = Dict(zip(varNames,[C1,C2,B1,B2,K1,K2,Ra,Rb]))
    return solDict, γ1, γ2
end

function NashRobustNation1gamma1(SOp,Dp,B1,B2,Rb,γ2)

    g(x) = (SOp.g∞ * x + SOp.g0) / (x+1)

    RobustProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(RobustProblem, "tol", GLOBAL_TOL)
    set_silent(RobustProblem)

    JuMP.@variables(RobustProblem, begin
    γ1o ≥ 0
    end)
   
    @NLexpression(RobustProblem, gRb, (SOp.g∞ * Rb + SOp.g0) / (Rb + 1))
    @NLexpression(RobustProblem, K2, 1/(SOp.A̅*SOp.h0-1) * ((SOp.A̅*SOp.h0-1) / (γ2*SOp.Φ*SOp.ηK) + B2) )

    @NLobjective(RobustProblem, Min,
    γ1o * (Dp.S̅/SOp.ρ - Dp.P0/SOp.ρ - Dp.T0/(SOp.ρ+SOp.ϕ) +
    - SOp.ηK*SOp.Φ/SOp.ρ * ((B1+Rb)/(SOp.A̅-1)+K2) + SOp.Φ*SOp.ηB/SOp.ρ * (B1^SOp.θ1 + gRb*B2^SOp.θ2)) +
    SOp.αR/SOp.ρ*(γ1o-SOp.γ̂1)^2 + 1/SOp.ρ*log((SOp.A̅-1)/(γ1o*SOp.Φ*SOp.ηK)) - 1/SOp.ρ)

    JuMP.optimize!(RobustProblem)

    γ1 = value(γ1o)
    
    return γ1
end

function NashRobusNation2gamma2(SOp,Dp,B1,B2,Rb,γ1)

    g(x) = (SOp.g∞ * x + SOp.g0) / (x+1)

    RobustProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(RobustProblem, "tol", GLOBAL_TOL)
    set_silent(RobustProblem)

    JuMP.@variables(RobustProblem, begin
    γ2o ≥ 0
    end)
   
    @NLexpression(RobustProblem, gRb, (SOp.g∞ * Rb + SOp.g0) / (Rb + 1))
    @NLexpression(RobustProblem, K1, 1/(SOp.A̅-1) * ((SOp.A̅-1) / (γ1*SOp.Φ*SOp.ηK) + B1 + Rb))

    @NLobjective(RobustProblem, Min,
    γ2o * (Dp.S̅/SOp.ρ - Dp.P0/SOp.ρ - Dp.T0/(SOp.ρ+SOp.ϕ) +
    - SOp.ηK*SOp.Φ/SOp.ρ * (K1 + B2/(SOp.A̅*SOp.h0-1)) + SOp.Φ*SOp.ηB/SOp.ρ * (B1^SOp.θ1 + gRb*B2^SOp.θ2)) +
    SOp.αR/SOp.ρ*(γ2o-SOp.γ̂2)^2 + 1/SOp.ρ*log((SOp.A̅*SOp.h0-1)/(γ2o*SOp.Φ*SOp.ηK)) - 1/SOp.ρ)

    JuMP.optimize!(RobustProblem)

    γ2 = value(γ2o)
    
    return γ2
end


function optimNashExplicit(p; quiet = true)

    #check that we are in the case where there is an explicit solution 
    @assert p.α == 1

    function computeRb(p; quiet = quiet)
        RbProblem = Model(Ipopt.Optimizer)
        set_optimizer_attribute(RbProblem, "tol", 1e-16)
        JuMP.@variable(RbProblem, Rb ≥ 0)

        @NLexpression(RbProblem, gRb, (p.g∞ * Rb + p.g0) / (Rb + 1))
        @NLexpression(RbProblem, gPrimeRb, (p.g∞ - p.g0) / (Rb + 1)^2 )

        @NLobjective(RbProblem, Min, 
        (p.ηK / p.ηB - gRb^p.θ2 * gPrimeRb^(1-p.θ2) * p.θ2^p.θ2 * (p.A̅ - 1)^(1-p.θ2) * (p.A̅ * p.h0 - 1)^p.θ2)^2 )

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
    
    if p.ηK / p.ηB >= p.g0^p.θ2 * p.gPrime0^(1-p.θ2) * p.θ2^p.θ2 * (p.A̅ - 1)^(1-p.θ2) * (p.A̅ * p.h0 - 1)^p.θ2
        if quiet == false
            println("Explicit condition for Rb = 0 matched")
        end
        Rb = 0 
    else
        Rb = computeRb(p; quiet = quiet)
    end
    B2 = (((p.ηB * p.θ2 * g(Rb)) * (p.A̅ * p.h0 - 1)) / p.ηK)^(1/(1 - p.θ2))
    C1 = (p.γ1*p.Φ*p.ηK/(p.A̅-1))^(-1/p.σ1)
    B1 = (p.ηK/(p.ηB*p.θ1*(p.A̅-1)))^(1/(p.θ1-1))
    Ra = 0
    K1 = 1/(p.A̅ - 1) * (C1 + B1 + Rb)
    C2 = (p.γ2*p.Φ*p.ηK/(p.A̅*p.h0-1))^(-1/p.σ2)
    K2 = 1/(p.A̅*p.h0 - 1) * (C2 + B2)
    
    solValues = [C1,C2,B1,B2,K1,K2,Ra,Rb]

    solDict = Dict(zip(varNames,solValues))
    return solDict
 
end

function computeGNash(p)
    solDict = optimNash(p; quiet = true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [solDict[name] for name in varNames]
    return G(K1,K2,B1,B2,Rb,p)
end

function computeWelFare1Nash(p,sol)

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

function computeWelFare2Nash(p,sol)

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


function computeAllNash(SOp,Dp; RaImpact = 1.0, RbImpact = 1.0)

    sol = optimNash(SOp,quiet=true,RaImpact=RaImpact,RbImpact=RbImpact)
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

    WelFare1 = computeWelFare1Nash(SOp,sol)
    WelFare2 = computeWelFare2Nash(SOp,sol)
    WelFare = WelFare1 + WelFare2
    Y1 = SOp.A̅ * K1^SOp.α
    Y2 = SOp.A̅ * h(Ra,SOp) * K2^SOp.α
    gRb = g(Rb,SOp)
    hRa = h(Ra,SOp)

    return C1,C2,B1,B2,K1,K2,Ra,Rb,gRb,hRa, GComputed,G1,G2, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end


function computeAllNashExplicit(SOp,Dp)

    sol = optimNashExplicit(SOp,quiet=true)
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

    WelFare1 = computeWelFare1Nash(SOp,sol)
    WelFare2 = computeWelFare2Nash(SOp,sol)
    WelFare = WelFare1 + WelFare2
    Y1 = SOp.A̅ * K1^SOp.α
    Y2 = SOp.A̅ * h(Ra,SOp) * K2^SOp.α
    gRb = g(Rb,SOp)
    hRa = h(Ra,SOp)

    return C1,C2,B1,B2,K1,K2,Ra,Rb,gRb,hRa, GComputed,G1,G2, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end
