function optimStackelberg(p; quiet = false)
    
    function Br2(x...)

        B1v,I1v,Rav,Rbv = x
        Br2Problem = Model(Ipopt.Optimizer)
        set_optimizer_attribute(Br2Problem, "tol", 1e-16)
    
        JuMP.@variables(Br2Problem, begin
            C2 ≥ 0.0
            B2 ≥ 0.0
            I2 ≥ 0.0
        end)

        JuMP.@NLparameters(Br2Problem, begin
            B1 == B1v 
            I1 == I1v
            Ra == Rav
            Rb == Rbv
        end) 
    
        # subexpressions
        @NLexpressions(Br2Problem, begin
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
            @NLexpression(Br2Problem, u2C2, log(C2))
        else
            @NLexpression(Br2Problem, u2C2, C2^(1-p.σ2)/(1-p.σ2))
        end
    
        # budget constraint
        @NLconstraint(Br2Problem,BudgetConstraint,
        C2 + B2 + I2  ≤ p.A̅ * hRa * fI2)
    
        # objective
        @NLobjective(Br2Problem, Max, u2C2 -  p.γ2 * p.Φ * GID)
    
        set_silent(Br2Problem)
        optimize!(Br2Problem)
    
        solValues = value.([C2, B2, I2])
        solValues[solValues .< 0] .= 0
        C2, B2, I2 = solValues
    
        return [C2, B2, I2]
    end
    
    Br2C2(x...) = Br2(x...)[1]
    Br2B2(x...) = Br2(x...)[2]
    Br2I2(x...) = Br2(x...)[3]

    function ∇Br2C2(g,x...)
        g .= grad(central_fdm(8, 1; factor=1e10),Br2C2, x...)
        return
    end

    function ∇Br2B2(g,x...)
        g .= grad(central_fdm(8, 1; factor=1e10),Br2B2, x...)
        return
    end

    function ∇Br2I2(g,x...)
        g .= grad(central_fdm(8, 1; factor=1e10),Br2I2, x...)
        return
    end

    stackelbergProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(stackelbergProblem, "tol", 1e-16)
    
    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    JuMP.@variables(stackelbergProblem,begin
        C1 ≥ 0.0
        B1 ≥ 0.0
        I1 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    register(stackelbergProblem, :Br2C2, 4, Br2C2, ∇Br2C2)
    register(stackelbergProblem, :Br2B2, 4, Br2B2, ∇Br2B2)
    register(stackelbergProblem, :Br2I2, 4, Br2I2, ∇Br2I2)

    @NLexpressions(stackelbergProblem, begin
        gRb, (p.g∞ * Rb + p.g0) / (Rb + 1)
        hRa, (p.h∞ * Ra + p.h0) / (Ra + 1)
        fI1, I1^p.α
        fI2, (Br2I2(B1,I1,Ra,Rb))^p.α
        D1, p.b̅ * B1^p.θ1
        D2, p.b̅ * gRb *  (Br2B2(B1,I1,Ra,Rb))^p.θ2
        I, I1 + Br2I2(B1,I1,Ra,Rb)
        D, D1 + D2
        GID, p.cI * I^p.δI - p.cD * D^p.δD
    end)

    if p.σ1 == 1
        @NLexpression(stackelbergProblem, u1C1, log(C1))
    else
        @NLexpression(stackelbergProblem, u1C1, C1^(1-p.σ1)/(1-p.σ1))
    end

    @NLconstraint(stackelbergProblem,BudgetConstraint1,
    C1 + B1 + I1 + Ra + Rb ≤ p.A̅ * fI1 )

    @NLobjective(stackelbergProblem, Max, u1C1 - p.γ1 * p.Φ * GID) 
   
    # verbose
    if quiet == false
        unset_silent(stackelbergProblem)
    else
        set_silent(stackelbergProblem)
    end

    optimize!(stackelbergProblem)

    C1v,B1v,I1v,Rav,Rbv = value.([C1,B1,I1,Ra,Rb])

    C2v,B2v,I2v = Br2(B1v,I1v,Rav,Rbv)
    solValues = [C1v,C2v,B1v,B2v,I1v,I2v,Rav,Rbv]    
    solValues[solValues .< 0] .= 0

    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    solDict = Dict(zip(varNames,solValues))
    return solDict
end
    

function optimStackelbergExplicit(p; quiet = true)

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
            # ricontrollare
            @NLobjective(RaRbProblem, Max, 
            - p.cI/(p.A̅ - 1) * (Ra + Rb) - (p.cI/(p.A̅ * hRa - 1))^(1-1/p.σ2) * (p.γ2 * p.Φ)^(-1/p.σ2) + 
            (p.cI/(p.A̅ * hRa - 1))^(-p.θ2/(1-p.θ2)) * (p.cD * gRb)^(1/(1-p.θ2))*p.θ2^(p.θ2/(1-p.θ2))*(1-p.θ2) )
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
    if p.cI / p.cD >= p.g0^p.θ2 * p.gPrime0^(1-p.θ2) * p.θ2^p.θ2 * (p.A̅ - 1)^(1-p.θ2) * (p.A̅ * h(Ra) - 1)^p.θ2
        println("Explicit condition for Rb = 0 matched")
    end

    C1 = ( (p.A̅ - 1) / (p.γ1*p.Φ*p.cI) )^(1/p.σ1)
    C2 = ((p.A̅ * h(Ra) - 1)/(p.γ2 * p.Φ * p.cI))^(1/p.σ2)
    B1 = (p.cD / p.cI * p.θ1 * (p.A̅ - 1))^(1/(1-p.θ1))
    B2 = (p.cD * g(Rb) * p.θ2 / p.cI * (p.A̅ * h(Ra) - 1))^(1/(1 - p.θ2))
    I1 = 1/(p.A̅ - 1) * (C1 + B1 + Ra + Rb)
    I2 = 1/(p.A̅*h(Ra) - 1) * (C2 + B2)

    
    solValues = [C1,C2,B1,B2,I1,I2,Ra,Rb]
    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    solDict = Dict(zip(varNames,solValues))
    
    return solDict
 
end

