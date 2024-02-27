function optimNash(p; TOL = 1e-6, MAX_IT = 1e6, quiet = true)
    # by iterated best response

    
    # initial point by social planner optimum
    solPlanner = optimPlanner(p; quiet = true, showObj = false)
    B2,K2 = solPlanner["B2"],solPlanner["K2"]
    
    C1,B1,K1,Ra,Rb = iterate1(p,B2,K2)
    C2, B2, K2 = iterate2(p,B1,K1,Ra,Rb)
    candidateSol = [C1,C2,B1,B2,K1,K2,Ra,Rb]

    err = 1.0; Nit = 0
    while (err > TOL) & (Nit < MAX_IT)
        
        C1,B1,K1,Ra,Rb = iterate1(p,B2,K2)
        C2, B2, K2 = iterate2(p,B1,K1,Ra,Rb)
        
        err = norm(candidateSol - [C1,C2,B1,B2,K1,K2,Ra,Rb],Inf)
        Nit = Nit + 1 

        candidateSol = [C1,C2,B1,B2,K1,K2,Ra,Rb]
        
        if quiet == false @show err, Nit end
    end
    
    if Nit == MAX_IT
        printstyled("nashProblem -> maximum iteration reached" * "\n",color=:red)
    end

    solDict = Dict(zip(varNames,candidateSol))

    return solDict
end

function iterate1(p,B2,K2)

    nash1Problem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(nash1Problem, "tol", 1e-32)
    set_optimizer_attribute(nash1Problem, "acceptable_tol",  1e-32)
    set_optimizer_attribute(nash1Problem, "max_iter", Int(1e6))

    JuMP.@variables(nash1Problem,begin
        C1 ≥ 0.0
        B1 ≥ 0.0
        K1 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(nash1Problem, begin
    gRb, (p.g∞ * Rb + p.g0) / (Rb + 1)
    hRa, (p.h∞ * Ra + p.h0) / (Ra + 1)
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
    solValues[solValues .< 0] .= 0
    C1,B1,K1,Ra,Rb = solValues

    return C1,B1,K1,Ra,Rb
end

function iterate2(p,B1,K1,Ra,Rb)

    nash2Problem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(nash2Problem, "tol", 1e-32)
    set_optimizer_attribute(nash2Problem, "acceptable_tol",  1e-32)
    set_optimizer_attribute(nash2Problem, "max_iter", Int(1e6))

    JuMP.@variables(nash2Problem,begin
        C2 ≥ 0.0
        B2 ≥ 0.0
        K2 ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(nash2Problem, begin
    gRb, (p.g∞ * Rb + p.g0) / (Rb + 1)
    hRa, (p.h∞ * Ra + p.h0) / (Ra + 1)
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
    solValues[solValues .< 0] .= 0
    C2, B2, K2 = solValues

    return C2, B2, K2
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
        println("Explicit condition for Rb = 0 matched")
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


function computeAllNash()

    sol = optimNash(SOpG,quiet=true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    GComputed = computeG(SOpG,sol)


    solODE = solveODE(SOpG,DpG,GComputed)
    P,T = solODE.u[end]
    TempFinal = ComputeTemperature(DpG,P,T)
    TempInitial = ComputeTemperature(DpG,DpG.P0,DpG.T0)

    ΔP = P - DpG.P0
    ΔT = T - DpG.T0
    ΔTemp = TempFinal - TempInitial

    WelFare1 = computeWelFare1Nash(SOpG,sol)
    WelFare2 = computeWelFare2Nash(SOpG,sol)
    WelFare = WelFare1 + WelFare2
    Y1 = SOpG.A̅ * K1^SOpG.α
    Y2 = SOpG.A̅ * h(Ra,SOpG) * K2^SOpG.α

    return C1,C2,B1,B2,K1,K2,Ra,Rb, GComputed, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end


function computeAllNashExplicit()

    sol = optimNashExplicit(SOpG,quiet=true)
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

    WelFare1 = computeWelFare1Nash(SOpG,sol)
    WelFare2 = computeWelFare2Nash(SOpG,sol)
    WelFare = WelFare1 + WelFare2
    Y1 = SOpG.A̅ * K1^SOpG.α
    Y2 = SOpG.A̅ * h(Ra,SOpG) * K2^SOpG.α
    gRb = g(Rb,SOpG)
    hRa = h(Ra,SOpG)

    return C1,C2,B1,B2,K1,K2,Ra,Rb,gRb,hRa, GComputed,G1,G2, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end
