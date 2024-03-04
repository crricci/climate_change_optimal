function optimNashRobust(SOp,Dp; quiet = false)
    SOpInner = deepcopy(SOp)

    #check that we are in the case where there is a soluton, i.e. AK log case
    @assert SOp.α == 1
    @assert SOp.σ1 == 1
    @assert SOp.σ2 == 1

    # initialization
    γ1 = SOpInner.γ1
    γ2 = SOpInner.γ2
    solNashStart = optimNashExplicit(SOpInner; quiet = true)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [solNashStart[name] for name in varNames]

    err = 1.0; Nit = 0
    while (err > GLOBAL_TOL) & (Nit < GLOBAL_MAX_IT)
        
        C1Pre,C2Pre,B1Pre,B2Pre,K1Pre,K2Pre,RaPre,RbPre,γ1Pre,γ2Pre = C1,C2,B1,B2,K1,K2,Ra,Rb,γ1,γ2 

        C1,B1,K1,Ra,Rb,γ1 = NashRobust1(SOpInner,Dp,C2,B2,K2,γ2; quiet = quiet)
        C2,B2,K2,γ2 = NashRobust2(SOpInner,Dp,C1,B1,K1,Ra,Rb,γ1; quiet = quiet)
        
        err = norm([C1,C2,B1,B2,K1,K2,Ra,Rb,γ1,γ2] - [C1Pre,C2Pre,B1Pre,B2Pre,K1Pre,K2Pre,RaPre,RbPre,γ1Pre,γ2Pre],Inf)
        Nit = Nit + 1 

        if quiet == false 
            print("Nash Case: ")
            @show err, Nit
         end
    end
    
    if Nit == GLOBAL_MAX_IT
        printstyled("nashProblem -> maximum iteration reached" * "\n",color=:red)
    end

    solDict = Dict(zip(varNames,[C1,C2,B1,B2,K1,K2,Ra,Rb]))
    return solDict, γ1, γ2
end

function NashRobust1(SOp,Dp,C2,B2,K2,γ2; quiet = false)

    SOpInner = deepcopy(SOp)
    SOpInner.γ2 = γ2

    # initial point by social planner optimum
    solNashStart = optimNashExplicit(SOpInner; quiet = true)
    C1,B1,K1,Ra,Rb = solNashStart["C1"],solNashStart["B1"],solNashStart["K1"],solNashStart["Ra"],solNashStart["Rb"]
    γ1 = 0.0

    err = 1.0; Nit = 0
    while (err > GLOBAL_TOL) & (Nit < GLOBAL_MAX_IT)
        
        γ1Pre = SOpInner.γ1
        C1Pre, B1Pre, K1Pre, RaPre, RbPre = C1,B1,K1,Ra,Rb

        C1,B1,K1,Ra,Rb = NashNation1(SOpInner,B2,K2)
        γ1 = NashRobust1Nature1(SOpInner,Dp,B1,B2,K2,Rb)
        
        err = norm([C1,B1,K1,Ra,Rb,γ1] - [C1Pre,B1Pre,K1Pre,RaPre,RbPre,γ1Pre],Inf)
        Nit = Nit + 1 

        SOpInner.γ1 = γ1

        if quiet == false 
            print("Nash Case Inner 1: ")
            @show err, Nit
         end        
    end
    
    if Nit == GLOBAL_MAX_IT
        printstyled("nashProblem -> maximum iteration reached" * "\n",color=:red)
    end

    return C1, B1, K1, Ra, Rb, γ1
end

function NashRobust1Nature1(SOp,Dp,B1,B2,K2,Rb)
    
    g(x) = (SOp.g∞ * x + SOp.g0) / (x+1)

    RobustProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(RobustProblem, "tol", GLOBAL_TOL)
    set_silent(RobustProblem)

    JuMP.@variables(RobustProblem, begin
    γ1o ≥ 0
    end)
   
    @NLexpression(RobustProblem, gRb, (SOp.g∞ * Rb + SOp.g0) / (Rb + 1))

    @NLobjective(RobustProblem, Min,
    γ1o * (Dp.S̅/SOp.ρ - Dp.P0/SOp.ρ - Dp.T0/(SOp.ρ+SOp.ϕ) +
    - SOp.ηK*SOp.Φ/SOp.ρ * ((B1+Rb)/(SOp.A̅-1)+K2) + SOp.Φ*SOp.ηB/SOp.ρ * (B1^SOp.θ1 + gRb*B2^SOp.θ2)) +
    SOp.αR/SOp.ρ*(γ1o-SOp.γ̂1)^2 + 1/SOp.ρ*log((SOp.A̅-1)/(γ1o*SOp.Φ*SOp.ηK)) - 1/SOp.ρ)

    JuMP.optimize!(RobustProblem)

    γ1 = value(γ1o)
    
    return γ1
end

function NashRobust2(SOp,Dp,C1,B1,K1,Ra,Rb,γ1; quiet = false)
    SOpInner = deepcopy(SOp)
    SOpInner.γ1 = γ1

    # initial point by social planner optimum
    solNashStart = optimNashExplicit(SOpInner; quiet = true)
    C2,B2,K2 = solNashStart["C2"],solNashStart["B2"],solNashStart["K2"]
    γ2 = 0.0

    err = 1.0; Nit = 0
    while (err > GLOBAL_TOL) & (Nit < GLOBAL_MAX_IT)
        
        γ2Pre = SOpInner.γ2
        C2Pre, B2Pre, K2Pre = C2, B2, K2

        C2, B2, K2 = NashNation2(SOpInner,B1,K1,Ra,Rb)
        γ2 = NashRobust1Nature2(SOpInner,Dp,B1,B2,K1,Rb)
        
        err = norm([C2,B2,K2,γ2] - [C2Pre,B2Pre,K2Pre,γ2Pre],Inf)
        Nit = Nit + 1 

        SOpInner.γ2 = γ2

        if quiet == false 
            print("Nash Case Inner 2: ")
            @show err, Nit
         end      
    end
    
    if Nit == GLOBAL_MAX_IT
        printstyled("nashProblem -> maximum iteration reached" * "\n",color=:red)
    end

    return C2,B2,K2,γ2
end

function NashRobust1Nature2(SOp,Dp,B1,B2,K1,Rb)
    g(x) = (SOp.g∞ * x + SOp.g0) / (x+1)

    RobustProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(RobustProblem, "tol", GLOBAL_TOL)
    set_silent(RobustProblem)

    JuMP.@variables(RobustProblem, begin
    γ2o ≥ 0
    end)
   
    @NLexpression(RobustProblem, gRb, (SOp.g∞ * Rb + SOp.g0) / (Rb + 1))

    @NLobjective(RobustProblem, Min,
    γ2o * (Dp.S̅/SOp.ρ - Dp.P0/SOp.ρ - Dp.T0/(SOp.ρ+SOp.ϕ) +
    - SOp.ηK*SOp.Φ/SOp.ρ * (K1 + B2/(SOp.A̅*SOp.h0-1)) + SOp.Φ*SOp.ηB/SOp.ρ * (B1^SOp.θ1 + gRb*B2^SOp.θ2)) +
    SOp.αR/SOp.ρ*(γ2o-SOp.γ̂2)^2 + 1/SOp.ρ*log((SOp.A̅*SOp.h0-1)/(γ2o*SOp.Φ*SOp.ηK)) - 1/SOp.ρ)

    JuMP.optimize!(RobustProblem)

    γ2 = value(γ2o)
    
    return γ2
end