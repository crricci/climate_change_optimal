
function optimPlannerNoResources(p; quiet = false)
    
    plannerNRProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerNRProblem, "tol", 1e-16)

    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    @variables(plannerNRProblem,begin
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
    @NLexpressions(plannerNRProblem, begin
    gRb, (Rb + p.g0) / (p.g∞ * Rb + 1)
    hRa, (Ra + p.h0) / (p.h∞ * Ra + 1)
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
    C1 + B1 + I1 + Ra + Rb ≤ p.A̅ * fI1)

    # budget constraint 2
    @NLconstraint(plannerNRProblem,BudgetConstraint2,
    C2 + B2 + I2 ≤ p.A̅ * hRa * fI2)
    
    Φ = (p.ϕL / p.ρ + (1-p.ϕL) * p.ϕ0 / (p.ρ + p.ϕ))

    # objective
    @NLobjective(plannerNRProblem, Max, u1C1 + u2C2 - (p.γ1 + p.γ2) * Φ * GID)

    # verbose
    if quiet == false
        unset_silent(plannerNRProblem)
    else
        set_silent(plannerNRProblem)
    end

    optimize!(plannerNRProblem)
    
    solValues = value.([C1,C2,B1,B2,I1,I2,Ra,Rb])
    solValues[solValues .< 0] .= 0
    solDict = Dict(zip(varNames,solValues))
    return solDict
end

