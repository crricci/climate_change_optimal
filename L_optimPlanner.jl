

function optimPlanner(p; quiet = false)
    
    plannerProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerProblem, "tol", 1e-16)

    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    JuMP.@variables(plannerProblem,begin
        C1 ≥ 0.0
        C2 ≥ 0.0
        B1 ≥ 0.0
        B2 ≥ 0.0
        I1 ≥ 0.0
        I2 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(plannerProblem, begin
    gRb, (p.g∞ * Rb + p.g0) / (Rb + 1)
    hRa, (p.h∞ * Ra + p.h0) / (Ra + 1)
    fI1, I1^p.α
    fI2, I2^p.α
    D1, p.b̅ * B1^p.θ1
    D2, p.b̅ * gRb * B2^p.θ2
    I, I1 + I2
    D, D1 + D2
    GID, p.cI * I^p.δI - p.cD * D^p.δD
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
    C1 + C2 + B1 + B2 + I1 + I2 + Ra + Rb ≤ p.A̅ * (fI1 + hRa * fI2))

    Φ = (p.ϕL / p.ρ + (1-p.ϕL) * p.ϕ0 / (p.ρ + p.ϕ))

    # objective
    @NLobjective(plannerProblem, Max, u1C1 + u2C2 - (p.γ1 + p.γ2) * Φ * GID)

    # verbose
    if quiet == false
        unset_silent(plannerProblem)
    else
        set_silent(plannerProblem)
    end

    optimize!(plannerProblem)
    
    solValues = value.([C1,C2,B1,B2,I1,I2,Ra,Rb])
    solValues[solValues .< 0] .= 0
    solDict = Dict(zip(varNames,solValues))
    return solDict
end

