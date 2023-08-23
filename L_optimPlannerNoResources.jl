
function optimPlannerNoResources(p; quiet = true, start = nothing)
    
    plannerNRProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerNRProblem, "tol", 1e-16)

    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    JuMP.@variables(plannerNRProblem,begin
        C1 ≥ 0.0
        C2 ≥ 0.0
        B1 ≥ 0.0
        B2 ≥ 0.0
        I1 ≥ 0.0
        I2 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    if !isnothing(start)
        set_start_value(C1,start["C1"])
        set_start_value(C2,start["C2"])
        set_start_value(B1,start["B1"])
        set_start_value(B2,start["B2"])
        set_start_value(I1,start["I1"])
        set_start_value(I2,start["I2"])
        set_start_value(Ra,start["Ra"])
        set_start_value(Rb,start["Rb"])
    end

    # subexpressions
    @NLexpressions(plannerNRProblem, begin
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
    
    # objective
    @NLobjective(plannerNRProblem, Max, u1C1 + u2C2 - (p.γ1 + p.γ2) * p.Φ * GID)

    # verbose
    if quiet == false
        unset_silent(plannerNRProblem)
    else
        set_silent(plannerNRProblem)
    end

    optimize!(plannerNRProblem)
    println("Objective: ",objective_value(plannerNRProblem))

    solValues = value.([C1,C2,B1,B2,I1,I2,Ra,Rb])
    solValues[solValues .< 0] .= 0
    solDict = Dict(zip(varNames,solValues))
    return solDict
end

function optimPlannerNoResourcesExplicit(p; quiet = true)

    #check that we are in the case where there is an explicit solution 
    @assert p.α == 1
    @assert p.δD == 1
    @assert p.δI == 1
    @assert p.b̅ == 1

    function computeRaRb(p; quiet = quiet)
        RaRbProblem = Model(Ipopt.Optimizer)
        set_optimizer_attribute(RaRbProblem, "tol", 1e-16)
        JuMP.@variable(RaRbProblem, Ra ≥ 0)
        JuMP.@variable(RaRbProblem, Rb ≥ 0)


        @NLexpression(RaRbProblem, gRb, (p.g∞ * Rb + p.g0) / (Rb + 1))
        @NLexpression(RaRbProblem, hRa, (p.h∞ * Ra + p.h0) / (Ra + 1))
        if p.σ2 == 1
            @NLobjective(RaRbProblem, Max, 
            1/((p.γ1+p.γ2)*p.Φ)*log((p.A̅*hRa - 1)/((p.γ1+p.γ2)*p.Φ*p.cI)) + 
            (1-p.θ2)*p.cD*gRb * (p.cI/ ((p.cD * p.θ2 * gRb) * (p.A̅*hRa - 1)))^(p.θ2/(p.θ2 - 1)) - p.cI/(p.A̅-1)*(Ra+Rb) )
        else
            @NLobjective(RaRbProblem, Max, 
            p.σ2/(1-p.σ2) * (1/((p.γ1+p.γ2)*p.Φ*p.cI))^(1/p.σ2) * (p.A̅*hRa - 1)^((1-p.σ2)/p.σ2) + 
            (1-p.θ2)*p.cD*gRb * (p.cI/ ((p.cD * p.θ2 * gRb) * (p.A̅*hRa - 1)))^(p.θ2/(p.θ2 - 1)) -
            p.cI/(p.A̅-1)*(Ra+Rb) )
        end

        # verbose
        if quiet == false
            unset_silent(RaRbProblem)
        else
            set_silent(RaRbProblem)
        end
        optimize!(RaRbProblem)
        Rb = value(Rb); Rb = Rb < 0 ? 0 : Rb
        Ra = value(Ra); Ra = Ra < 0 ? 0 : Ra

        return Ra,Rb
    end

    g(x) = (p.g∞ * x + p.g0) / (x+1)
    h(x) = (p.h∞ * x + p.h0) / (x+1)
    
    Ra,Rb = computeRaRb(p; quiet=quiet)
    C1 = ( (p.γ1+p.γ2) * p.Φ * p.cI / (p.A̅ - 1) )^(-1/p.σ1)
    C2 = ((p.A̅ * h(Ra) - 1)/((p.γ1+p.γ2) * p.Φ * p.cI))^(1/p.σ2)
    B1 = (p.cD/p.cI * p.θ1 * (p.A̅ - 1))^(1/(1-p.θ1))
    B2 = (p.cI/ ((p.cD * p.θ2 * g(Rb)) * (p.A̅ * h(Ra) - 1)))^(1/(p.θ2 - 1))
    I1 = 1/(p.A̅ - 1) * (C1 + B1 + Ra + Rb)
    I2 = 1/(p.A̅*h(Ra) - 1) * (C2 + B2)

    
    solValues = [C1,C2,B1,B2,I1,I2,Ra,Rb]
    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    solDict = Dict(zip(varNames,solValues))

    println("Objective: ", objectiveValuePlannerNoResources(p,solDict))
    
    return solDict
 
end

function objectiveValuePlannerNoResources(p,sol)

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