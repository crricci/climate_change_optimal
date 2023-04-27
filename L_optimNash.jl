function optimNash(p; TOL = 1e-16, MAX_IT = 1e6, quiet = false)
    # by iterated best response

    
    # initial point by social planner optimum
    solPlanner = optimPlanner(p; quiet = true)
    B2,I2 = solPlanner["B2"],solPlanner["I2"]
    
    C1,B1,I1,Ra,Rb = iterate1(p,B2,I2)
    C2, B2, I2 = interate2(p,B1,I1,Ra,Rb)
    candidateSol = [C1,C2,B1,B2,I1,I2,Ra,Rb]

    err = 1.0; Nit = 0
    while (err > TOL) & (Nit < MAX_IT)
        
        C1,B1,I1,Ra,Rb = iterate1(p,B2,I2)
        C2, B2, I2 = interate2(p,B1,I1,Ra,Rb)
        
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

    @variables(nash1Problem,begin
        C1 ≥ 0.0
        B1 ≥ 0.0
        I1 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(nash1Problem, begin
    gRb, (Rb + p.g0) / (p.g∞ * Rb + 1)
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

    Φ = (p.ϕL / p.ρ + (1-p.ϕL) * p.ϕ0 / (p.ρ + p.ϕ))

    # objective
    @NLobjective(nash1Problem, Max, u1C1 - p.γ1 * Φ * GID)

    set_silent(nash1Problem)
    optimize!(nash1Problem)
    
    solValues = value.([C1,B1,I1,Ra,Rb])
    solValues[solValues .< 0] .= 0
    C1,B1,I1,Ra,Rb = solValues

    return C1,B1,I1,Ra,Rb
end

function interate2(p,B1,I1,Ra,Rb)

    nash2Problem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(nash2Problem, "tol", 1e-16)

    @variables(nash2Problem,begin
        C2 ≥ 0.0
        B2 ≥ 0.0
        I2 ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(nash2Problem, begin
    gRb, (Rb + p.g0) / (p.g∞ * Rb + 1)
    hRa, (Ra + p.h0) / (p.h∞ * Ra + 1)
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

    Φ = (p.ϕL / p.ρ + (1-p.ϕL) * p.ϕ0 / (p.ρ + p.ϕ))

    # objective
    @NLobjective(nash2Problem, Max, u2C2 -  p.γ2 * Φ * GID)

    set_silent(nash2Problem)
    optimize!(nash2Problem)

    solValues = value.([C2, B2, I2])
    solValues[solValues .< 0] .= 0
    C2, B2, I2 = solValues

    return C2, B2, I2
end