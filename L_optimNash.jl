function optimNash(p; TOL = 1e-16, MAX_IT = 1e6, quiet = true)
    # by iterated best response

    
    # initial point by social planner optimum
    solPlanner = optimPlanner(p; quiet = true, showObj = false)
    B2,I2 = solPlanner["B2"],solPlanner["I2"]
    
    C1,B1,I1,Ra,Rb = iterate1(p,B2,I2)
    C2, B2, I2 = iterate2(p,B1,I1,Ra,Rb)
    candidateSol = [C1,C2,B1,B2,I1,I2,Ra,Rb]

    err = 1.0; Nit = 0
    while (err > TOL) & (Nit < MAX_IT)
        
        C1,B1,I1,Ra,Rb = iterate1(p,B2,I2)
        C2, B2, I2 = iterate2(p,B1,I1,Ra,Rb)
        
        err = norm(candidateSol - [C1,C2,B1,B2,I1,I2,Ra,Rb],Inf)
        Nit = Nit + 1 

        candidateSol = [C1,C2,B1,B2,I1,I2,Ra,Rb]
        
        if quiet == false @show err, Nit end
    end
    
    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    solDict = Dict(zip(varNames,candidateSol))

    return solDict
end

function iterate1(p,B2,I2)

    nash1Problem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(nash1Problem, "tol", 1e-16)

    JuMP.@variables(nash1Problem,begin
        C1 ≥ 0.0
        B1 ≥ 0.0
        I1 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(nash1Problem, begin
    gRb, (p.g∞ * Rb + p.g0) / (Rb + 1)
    hRa, (p.h∞ * Ra + p.h0) / (Ra + 1)
    fI1, I1^p.α
    D1, p.b̅ * B1^p.θ1
    D2, p.b̅ * gRb * B2^p.θ2
    I, I1 + I2
    D, D1 + D2
    GID, p.cI * I^p.δI - p.cD * D^p.δD
    end)
    
    # utility
    if p.σ1 == 1
        @NLexpression(nash1Problem, u1C1, log(C1))
    else
        @NLexpression(nash1Problem, u1C1, C1^(1-p.σ1)/(1-p.σ1))
    end

    # budget constraint
    @NLconstraint(nash1Problem,BudgetConstraint,
    C1 + B1 + I1 + Ra + Rb ≤ p.A̅ * fI1)

    # objective
    @NLobjective(nash1Problem, Max, u1C1 - p.γ1 * p.Φ * GID)

    set_silent(nash1Problem)
    optimize!(nash1Problem)
    
    solValues = value.([C1,B1,I1,Ra,Rb])
    solValues[solValues .< 0] .= 0
    C1,B1,I1,Ra,Rb = solValues

    return C1,B1,I1,Ra,Rb
end

function iterate2(p,B1,I1,Ra,Rb)

    nash2Problem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(nash2Problem, "tol", 1e-16)

    JuMP.@variables(nash2Problem,begin
        C2 ≥ 0.0
        B2 ≥ 0.0
        I2 ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(nash2Problem, begin
    gRb, (p.g∞ * Rb + p.g0) / (Rb + 1)
    hRa, (p.h∞ * Ra + p.h0) / (Ra + 1)
    fI2, I2^p.α
    D1, p.b̅ * B1^p.θ1
    D2, p.b̅ * gRb * B2^p.θ2
    I, I1 + I2
    D, D1 + D2
    GID, p.cI * I^p.δI - p.cD * D^p.δD
    end)
    
    # utility
    if p.σ2 == 1
        @NLexpression(nash2Problem, u2C2, log(C2))
    else
        @NLexpression(nash2Problem, u2C2, C2^(1-p.σ2)/(1-p.σ2))
    end

    # budget constraint
    @NLconstraint(nash2Problem,BudgetConstraint,
    C2 + B2 + I2  ≤ p.A̅ * hRa * fI2)

    # objective
    @NLobjective(nash2Problem, Max, u2C2 -  p.γ2 * p.Φ * GID)

    set_silent(nash2Problem)
    optimize!(nash2Problem)

    solValues = value.([C2, B2, I2])
    solValues[solValues .< 0] .= 0
    C2, B2, I2 = solValues

    return C2, B2, I2
end

function optimNashExplicit(p; quiet = true)

    #check that we are in the case where there is an explicit solution 
    @assert p.α == 1
    @assert p.δD == 1
    @assert p.δI == 1
    @assert p.b̅ == 1

    function computeRb(p; quiet = quiet)
        RbProblem = Model(Ipopt.Optimizer)
        set_optimizer_attribute(RbProblem, "tol", 1e-16)
        JuMP.@variable(RbProblem, Rb ≥ 0)

        @NLexpression(RbProblem, gRb, (p.g∞ * Rb + p.g0) / (Rb + 1))
        @NLexpression(RbProblem, gPrimeRb, (p.g∞ - p.g0) / (Rb + 1)^2 )

        @NLobjective(RbProblem, Min, 
        (p.cI / p.cD - gRb^p.θ2 * gPrimeRb^(1-p.θ2) * p.θ2^p.θ2 * (p.A̅ - 1)^(1-p.θ2) * (p.A̅ * p.h0 - 1)^p.θ2)^2 )

        # verbose
        if quiet == false
            unset_silent(RbProblem)
        else
            set_silent(RbProblem)
        end

        optimize!(RbProblem)
        Rb = value(Rb); Rb = Rb < 0 ? 0 : Rb

        return Rb
    end

    g(x) = (p.g∞ * x + p.g0) / (x+1)
    
    if p.cI / p.cD >= p.g0^p.θ2 * p.gPrime0^(1-p.θ2) * p.θ2^p.θ2 * (p.A̅ - 1)^(1-p.θ2) * (p.A̅ * p.h0 - 1)^p.θ2
        println("Explicit condition for Rb = 0 matched")
        Rb = 0 
    else
        Rb = computeRb(p; quiet = quiet)
    end
    B2 = (((p.cD * p.θ2 * g(Rb)) * (p.A̅ * p.h0 - 1)) / p.cI)^(1/(1 - p.θ2))
    C1 = (p.γ1*p.Φ*p.cI/(p.A̅-1))^(-1/p.σ1)
    B1 = (p.cI/(p.cD*p.θ1*(p.A̅-1)))^(1/(p.θ1-1))
    Ra = 0
    I1 = 1/(p.A̅ - 1) * (C1 + B1 + Rb)
    C2 = (p.γ2*p.Φ*p.cI/(p.A̅*p.h0-1))^(-1/p.σ2)
    I2 = 1/(p.A̅*p.h0 - 1) * (C2 + B2)
    
    solValues = [C1,C2,B1,B2,I1,I2,Ra,Rb]
    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]

    solDict = Dict(zip(varNames,solValues))
    return solDict
 
end

