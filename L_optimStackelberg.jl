
function optimStackelberg(p; quiet = true, showObj = false)

    stackelbergProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(stackelbergProblem, "tol", GLOBAL_TOL)
    set_optimizer_attribute(stackelbergProblem, "acceptable_tol",  GLOBAL_TOL)
    set_optimizer_attribute(stackelbergProblem, "max_iter", GLOBAL_MAX_IT)

    JuMP.@variables(stackelbergProblem, begin
    C1 ≥ 0.0
    B1 ≥ 0.0
    K1 ≥ 0.0
    Ra ≥ 0.0
    Rb ≥ 0.0

    # inner problem
    C2 ≥ 0.0    # primal feasibility 
    B2 ≥ 0.0    # primal feasibility 
    K2 ≥ 0.0    # primal feasibility 
    μC2 ≥ 0     # dual feasibility 
    μB2 ≥ 0     # dual feasibility 
    μK2 ≥ 0     # dual feasibility 
    λ2 ≥ 0      # dual feasibility 
    end)

    # subexpressions
    @NLexpressions(stackelbergProblem, begin
    gRb, (p.g∞ * Rb + p.g0) / (Rb + 1)
    hRa, (p.h∞ * Ra + p.h0) / (Ra + 1)
    fK1, K1^p.α
    fK2, K2^p.α
    D1, B1^p.θ1
    D2, gRb * B2^p.θ2
    K, K1 + K2
    D, D1 + D2
    GID, p.ηK * K - p.ηB * D
    # derivatives for inner problem
    ∂u2partialC2, 1/C2^p.σ2
    ∂G∂B2, -p.ηB * gRb * p.θ2 * B2^(p.θ2-1)
    ∂G∂K2, p.ηK
    ∂fK2∂K2, p.α * K2^(p.α-1)
    end)

    # budget constraint 1
    @NLconstraint(stackelbergProblem,BudgetConstraint1,
    C1 + B1 + K1 + Ra + Rb ≤ p.A̅ * fK1)

    # budget constraint 2, primal feasibility
    @NLconstraint(stackelbergProblem,BudgetConstraint2,
    C2 + B2 + K2 ≤ p.A̅ * hRa * fK2)

    # utility
    if p.σ1 == 1
        @NLexpression(stackelbergProblem, u1C1, log(C1))
    else
        @NLexpression(stackelbergProblem, u1C1, C1^(1-p.σ1)/(1-p.σ1))
    end
    if p.σ2 == 1
        @NLexpression(stackelbergProblem, u2C2, log(C2))
    else
        @NLexpression(stackelbergProblem, u2C2, C2^(1-p.σ2)/(1-p.σ2))
    end

    # stationarity
    @NLconstraints(stackelbergProblem, begin
    - ∂u2partialC2 + (-μC2 + λ2) == 0
    p.γ2 * p.Φ * ∂G∂B2 + (-μB2 + λ2) == 0
    p.γ2 * p.Φ * ∂G∂K2 + (-μK2 + λ2 -λ2 * p.A̅ * hRa * ∂fK2∂K2) == 0
    end)

    # complementary slackness
    # regularization
    t = 1e-20
    @NLconstraints(stackelbergProblem, begin
    μC2 * ( - C2) ≥ -t
    μB2 * ( - B2) ≥ -t
    μK2 * ( - K2) ≥ -t
    λ2 * (C2 + B2 + K2 - p.A̅ * hRa * fK2) ≥ -t
    end)

    # @NLconstraints(stackelbergProblem, begin
    # μC2 * ( - C2) == 0
    # μB2 * ( - B2) == 0
    # μK2 * ( - K2) == 0
    # λ2 * (C2 + B2 + K2 - p.A̅ * hRa * fK2) == 0
    # end)
    
    @NLobjective(stackelbergProblem, Max, u1C1 - p.γ1 * p.Φ * GID)

    if quiet == false
        unset_silent(stackelbergProblem)
    else
        set_silent(stackelbergProblem)
    end

    JuMP.optimize!(stackelbergProblem)

    if termination_status(stackelbergProblem) ∉ [LOCALLY_SOLVED, OPTIMAL]
        printstyled("stackelbergProblem, termination status -> " * string(termination_status(stackelbergProblem)) * "\n",color=:red)
    end

    if showObj println("Objective: ",objective_value(stackelbergProblem)) end

    solValues = value.([C1,C2,B1,B2,K1,K2,Ra,Rb])
    solValues[solValues .< 0] .= 0
    solDict = Dict(zip(varNames,solValues))
    return solDict
end
    

function optimStackelbergExplicit(p; quiet = true)

    #check that we are in the case where there is an explicit solution 
    @assert p.α == 1

    function computeRaRb(p; quiet = quiet)
        RaRbProblem = Model(Ipopt.Optimizer)
        # set_optimizer_attribute(RaRbProblem, "tol", 1e-16)
        JuMP.@variable(RaRbProblem, Ra ≥ 0)
        JuMP.@variable(RaRbProblem, Rb ≥ 0)


        @NLexpression(RaRbProblem, gRb, (p.g∞ * Rb + p.g0) / (Rb + 1))
        @NLexpression(RaRbProblem, hRa, (p.h∞ * Ra + p.h0) / (Ra + 1))

        # most likely the formula for sigma2 = 1 is wrong
        # the one for sigma2 != 1 is corect, seems also for sigma2 = 1
        
        # if p.σ2 == 1
            # @NLobjective(RaRbProblem, Max, 
            # 1/((p.γ1+p.γ2)*p.Φ)*log((p.A̅*hRa - 1)/((p.γ1+p.γ2)*p.Φ*p.ηK)) + 
            # (1-p.θ2)*p.ηB*gRb * (p.ηK/ ((p.ηB * p.θ2 * gRb) * (p.A̅*hRa - 1)))^(p.θ2/(p.θ2 - 1)) - p.ηK/(p.A̅-1)*(Ra+Rb) )
        # else
            @NLobjective(RaRbProblem, Max, 
            - p.ηK/(p.A̅ - 1) * (Ra + Rb) - (p.ηK/(p.A̅ * hRa - 1))^(1-1/p.σ2) * (p.γ2 * p.Φ)^(-1/p.σ2) + 
            (p.ηK/(p.A̅ * hRa - 1))^(-p.θ2/(1-p.θ2)) * (p.ηB * gRb)^(1/(1-p.θ2))*p.θ2^(p.θ2/(1-p.θ2))*(1-p.θ2) )
        # end

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
    if p.ηK / p.ηB >= p.g0^p.θ2 * p.gPrime0^(1-p.θ2) * p.θ2^p.θ2 * (p.A̅ - 1)^(1-p.θ2) * (p.A̅ * h(Ra) - 1)^p.θ2
        println("Explicit condition for Rb = 0 matched")
    end

    C1 = ( (p.A̅ - 1) / (p.γ1*p.Φ*p.ηK) )^(1/p.σ1)
    C2 = ((p.A̅ * h(Ra) - 1)/(p.γ2 * p.Φ * p.ηK))^(1/p.σ2)
    B1 = (p.ηB / p.ηK * p.θ1 * (p.A̅ - 1))^(1/(1-p.θ1))
    B2 = (p.ηB * g(Rb) * p.θ2 / p.ηK * (p.A̅ * h(Ra) - 1))^(1/(1 - p.θ2))
    K1 = 1/(p.A̅ - 1) * (C1 + B1 + Ra + Rb)
    K2 = 1/(p.A̅*h(Ra) - 1) * (C2 + B2)

    
    solValues = [C1,C2,B1,B2,K1,K2,Ra,Rb]
    solDict = Dict(zip(varNames,solValues))
    
    return solDict
 
end

function computeGStackelberg(p)

    solDict = optimStackelberg(p; quiet = true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [solDict[name] for name in varNames]
    return G(K1,K2,B1,B2,Rb,p)
end

function computeWelFare1Stackelberg(p,sol)

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

function computeWelFare2Stackelberg(p,sol)

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

function computeAllStackelberg()

    sol = optimStackelberg(SOpG,quiet=true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    GComputed = computeG(SOpG,sol)

    solODE = solveODE(SOpG,DpG,GComputed)
    P,T = solODE.u[end]
    TempFinal = ComputeTemperature(DpG,P,T)
    TempInitial = ComputeTemperature(DpG,DpG.P0,DpG.T0)

    ΔP = P - DpG.P0
    ΔT = T - DpG.T0
    ΔTemp = TempFinal - TempInitial

    WelFare1 = computeWelFare1Stackelberg(SOpG,sol)
    WelFare2 = computeWelFare2Stackelberg(SOpG,sol)
    WelFare = WelFare1 + WelFare2
    Y1 = SOpG.A̅ * K1^SOpG.α
    Y2 = SOpG.A̅ * h(Ra,SOpG) * K2^SOpG.α

    return C1,C2,B1,B2,K1,K2,Ra,Rb, GComputed, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end

function computeAllStackelbergExplicit()

    sol = optimStackelbergExplicit(SOpG,quiet=true)
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

    WelFare1 = computeWelFare1Stackelberg(SOpG,sol)
    WelFare2 = computeWelFare2Stackelberg(SOpG,sol)
    WelFare = WelFare1 + WelFare2
    Y1 = SOpG.A̅ * K1^SOpG.α
    Y2 = SOpG.A̅ * h(Ra,SOpG) * K2^SOpG.α
    gRb = g(Rb,SOpG)
    hRa = h(Ra,SOpG)

    return C1,C2,B1,B2,K1,K2,Ra,Rb,gRb,hRa, GComputed,G1,G2, ΔP, ΔT, ΔTemp, WelFare, WelFare1, WelFare2, Y1, Y2

end
