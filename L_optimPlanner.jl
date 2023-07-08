

function optimPlanner(p; quiet = true, showObj = true)
    
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

    # objective
    @NLobjective(plannerProblem, Max, u1C1 + u2C2 - (p.γ1 + p.γ2) * p.Φ * GID)

    # verbose
    if quiet == false
        unset_silent(plannerProblem)
    else
        set_silent(plannerProblem)
    end

    optimize!(plannerProblem)
    if showObj println("Objective: ",objective_value(plannerProblem)) end
    
    solValues = value.([C1,C2,B1,B2,I1,I2,Ra,Rb])
    solValues[solValues .< 0] .= 0
    solDict = Dict(zip(varNames,solValues))
    return solDict
end

function optimPlannerExplicit(p; quiet = false)

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
        @NLobjective(RbProblem, Max, -Rb + gRb^(1/(1-p.θ2)) * (p.θ2 / (p.cI / ((p.A̅ - 1) * p.cD)))^(1/(1-p.θ2)) * (1/p.θ2 - 1) )

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

    C1 = ( (p.A̅-1)/((p.γ1+p.γ2)*p.Φ*p.cI) )^(1/p.σ1)
    C2 = ( (p.A̅-1)/((p.γ1+p.γ2)*p.Φ*p.cI) )^(1/p.σ2)
    Ra = 0
    if p.cI / p.cD >= p.g0^p.θ2 * p.gPrime0^(1-p.θ2) * p.θ2^p.θ2 * (p.A̅ - 1)
        println("Explicit condition for Rb = 0 matched")
        Rb = 0 
    else
        Rb = computeRb(p; quiet = quiet)
    end
    B1 = (p.cD/p.cI * p.θ1 * (p.A̅ - 1))^(1/(1-p.θ1))
    B2 = (p.cD/p.cI * p.θ2 * g(Rb) * (p.A̅ - 1))^(1/(1-p.θ2))
    I1 = 1/(p.A̅ - 1) * (C1 + C2 + B1 + B2 + Rb)
    I2 = 0

    solDict = Dict("C1" => C1, "C2" => C2, "Ra" => Ra, "Rb" => Rb,
        "B1" => B1, "B2" => B2, "I1" => I1,"I2" => I2)

    return solDict

end

function objectiveValuePlanner(p,sol)

    B1, B2, C1, C2, I1, I2, Ra, Rb = values(sort(sol))

    # objective value
    u1C1 = p.σ1 == 1 ? log(C1) : C1^(1-p.σ1)/(1-p.σ1)
    u2C2 = p.σ2 == 1 ? log(C2) : C2^(1-p.σ2)/(1-p.σ2)
    gRb = (p.g∞ * Rb + p.g0) / (Rb + 1)
    hRa = (p.h∞ * Ra + p.h0) / (Ra + 1)
    D1 = p.b̅ * B1^p.θ1
    D2 = p.b̅ * gRb * B2^p.θ2
    I = I1 + I2
    D = D1 + D2
    GID = p.cI * I^p.δI - p.cD * D^p.δD

    obj = u1C1 + u2C2 - (p.γ1+p.γ2)*p.Φ*GID

    return obj
end