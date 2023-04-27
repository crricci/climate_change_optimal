

function optimPlanner(p)
    
    plannerProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerProblem, "tol", 1e-16)

    epsilon = 0.0
    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    @variables(plannerProblem,begin
        C1 ≥ epsilon
        C2 ≥ epsilon
        B1 ≥ 0 
        B2 ≥ 0 
        I1 ≥ 0
        I2 ≥ 0
        Ra ≥ 0
        Rb ≥ 0
    end)

    @NLexpressions(plannerProblem, begin
    gRb, (Rb + p.g0) / (p.g∞ * Rb + 1)
    hRa, (Ra + p.h0) / (p.h∞ * Ra + 1)
    fI1, I1^p.α
    fI2, I2^p.α
    D1, p.b̅ * B1^p.θ1
    D2, p.b̅ * gRb * B2^p.θ2

    # objective
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

    @NLconstraint(plannerProblem,BudgetConstraint,
    C1 + C2 + B1 + B2 + I1 + I2 + Ra + Rb ≤ p.A̅ * (fI1 + hRa * fI2))

    Φ = (p.γ1 + p.k * p.γ2) * (p.ϕL / p.ρ + (1-p.ϕL) * p.ϕ0 / (p.ρ + p.ϕ))

    @NLobjective(plannerProblem, Max, u1C1 + p.k * u2C2 - Φ*GID)

    unset_silent(plannerProblem)
    optimize!(plannerProblem)

    solution = value.([C1,C2,B1,B2,I1,I2,Ra,Rb])

    return(solution)
end

