

function FinalTempChangeSigma(SOp,Dp,NPtσ)


    # range of σ to test
    σRange = LinRange(0.1,1.9,NPtσ)

    TempFinalGlobalPlanner = zeros(NPtσ)
    TempFinalRestrictedPlanner = zeros(NPtσ)
    TempFinalNash = zeros(NPtσ)

    pTmp  = deepcopy(SOp)

    @showprogress for (i,σi) in enumerate(σRange)
        pTmp.σ2 = σi
        pTmp.σ1 = σi
        
        GlobalPlanner = optimPlannerExplicit(pTmp)
        C1,C2,B1,B2,K1,K2,Ra,Rb = [GlobalPlanner[name] for name in varNames]
        GComputed = computeG(pTmp,GlobalPlanner)
        G1 = computeG1(pTmp,GlobalPlanner)
        G2 = computeG2(pTmp,GlobalPlanner)
        solODE = solveODE(pTmp,Dp,GComputed)
        P,T = solODE.u[end]
        TempFinalGlobalPlanner[i] = ComputeTemperature(Dp,P,T)

        RestrictedPlanner = optimPlannerNoResourcesExplicit(pTmp)
        C1,C2,B1,B2,K1,K2,Ra,Rb = [RestrictedPlanner[name] for name in varNames]
        GComputed = computeG(pTmp,RestrictedPlanner)
        G1 = computeG1(pTmp,RestrictedPlanner)
        G2 = computeG2(pTmp,RestrictedPlanner)
        solODE = solveODE(pTmp,Dp,GComputed)
        P,T = solODE.u[end]
        TempFinalRestrictedPlanner[i] = ComputeTemperature(Dp,P,T)
        
        Nash = optimNashExplicit(pTmp)
        C1,C2,B1,B2,K1,K2,Ra,Rb = [Nash[name] for name in varNames]
        GComputed = computeG(pTmp,Nash)
        G1 = computeG1(pTmp,Nash)
        G2 = computeG2(pTmp,Nash)
        solODE = solveODE(pTmp,Dp,GComputed)
        P,T = solODE.u[end]
        TempFinalNash[i] = ComputeTemperature(Dp,P,T)
        
    end

    return σRange,TempFinalGlobalPlanner,TempFinalRestrictedPlanner,TempFinalNash
end

function plotFinalTempChangeSigma(SOp,Dp; NPtσ = 100)
    
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    sigmaRange,TempFinalGlobalPlanner,TempFinalRestrictedPlanner,TempFinalNash = FinalTempChangeSigma(SOp,Dp,NPtσ)

    plot(sigmaRange, TempFinalGlobalPlanner, label="Global Planner",color="blue")
    plot(sigmaRange, TempFinalRestrictedPlanner, label="Restricted Planner",color="gold")
    plot(sigmaRange, TempFinalNash, label="Nash",color="green")

    xlabel("σ")
    ylabel("Final Temperature (°C)")
    legend()
    grid(true)
    tight_layout()

    savefig("FinalTempChangeSigma.pdf")

end

function DeltaTempChangeSigmaRel(SOp,Dp,NPtσ)


    # range of σ to test
    σRange = LinRange(0.1,1.9,NPtσ)

    DeltaTemp = zeros(NPtσ)

    pTmp  = deepcopy(SOp)

    @showprogress for (i,σi) in enumerate(σRange)
        pTmp.σ2 = σi
        pTmp.σ1 = σi
        
        planner = optimPlannerExplicit(pTmp)
        C1,C2,B1,B2,K1,K2,Ra,Rb = [planner[name] for name in varNames]
        GComputed = computeG(pTmp,planner)
        G1 = computeG1(pTmp,planner)
        G2 = computeG2(pTmp,planner)
        solODE = solveODE(pTmp,Dp,GComputed)
        P,T = solODE.u[end]
        TempFinalPlanner = ComputeTemperature(Dp,P,T)
        
        nash = optimNashExplicit(pTmp)
        C1,C2,B1,B2,K1,K2,Ra,Rb = [nash[name] for name in varNames]
        GComputed = computeG(pTmp,nash)
        G1 = computeG1(pTmp,nash)
        G2 = computeG2(pTmp,nash)
        solODE = solveODE(pTmp,Dp,GComputed)
        P,T = solODE.u[end]
        TempFinalNash = ComputeTemperature(Dp,P,T)
        
        DeltaTemp[i] = (TempFinalNash-TempFinalPlanner)/TempFinalPlanner

    end

    return σRange,DeltaTemp
end

function plotDeltaTempChangeSigmaRel(SOp,Dp; NPtσ = 100)
    
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    sigmaRange, DeltaTemp = DeltaTempChangeSigma(SOp,Dp,NPtσ)

    plot(sigmaRange, 100 * DeltaTemp,color="black")
    plot([sigmaRange[findmax(DeltaTemp)[2]],sigmaRange[findmax(DeltaTemp)[2]]],100*[minimum(DeltaTemp),maximum(DeltaTemp)],)

    xlabel("σ")
    ylabel("Relative change in final temperature (%)")
    grid(true)
    tight_layout()

    savefig("DeltaTempChangeSigmaRel.pdf")

end

function RaRbChangeGamma(p,NPtγ2)


    # range of γ2 to test
    γ2Range = LinRange(p.γ2,2*p.γ2,NPtγ2)

    RaPlanner = zeros(NPtγ2)
    RbPlanner = zeros(NPtγ2)
    RaPlannerNoResources = zeros(NPtγ2)
    RbPlannerNoResources = zeros(NPtγ2)
    RaNash = zeros(NPtγ2)
    RbNash = zeros(NPtγ2)

    pTmp  = deepcopy(p)

    @showprogress for (i,γ2i) in enumerate(γ2Range)
        pTmp.γ2 = γ2i

        planner = optimPlannerExplicit(pTmp)
        RaPlanner[i] = planner["Ra"]
        RbPlanner[i] = planner["Rb"]

        plannerNoResources = optimPlannerNoResourcesExplicit(pTmp)
        RaPlannerNoResources[i] = plannerNoResources["Ra"]
        RbPlannerNoResources[i] = plannerNoResources["Rb"]

        nash = optimNashExplicit(pTmp)
        RaNash[i] = nash["Ra"]
        RbNash[i] = nash["Ra"]
    end

    return γ2Range,RaPlanner,RbPlanner,RaPlannerNoResources,RbPlannerNoResources,RaNash,RbNash

end

function plotRaRbChangeGamma(p; NPtγ2 = 10)
    
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    (γ2Range,RaPlanner,RbPlanner,RaPlannerNoResources,
        RbPlannerNoResources,RaNash,RbNash) = RaRbChangeGamma(p,NPtγ2)
    γ2Range = γ2Range/p.γ1


    # plot(γ2Range, RaPlanner/RaPlanner[1]*100, label = latexstring("\$R_a\$ Global Planner"),color="blue")
    plot(γ2Range, RaPlannerNoResources/RaPlannerNoResources[1]*100, label = latexstring("\$R_a\$ Restricted planner"),color="black")
    # plot(γ2Range, RaNash/RaNash[1]*100, label = latexstring("\$R_a\$ Nash"),color="green")

    # plot(γ2Range, RbPlanner/RbPlanner[1]*100, label = latexstring("\$R_b\$ Global Planner"),color="blue","--")
    plot(γ2Range, RbPlannerNoResources/RbPlannerNoResources[1]*100, label = latexstring("\$R_b\$ Restricted Planner"),color="black","--")
    # plot(γ2Range, RbNash/RbNash[1]*100, label = latexstring("\$R_b\$ Nash"),color="green","--")

    xlabel("γ2/γ1")
    ylabel(latexstring("Expenditure in Transfers (%) w.r.t. \$\\gamma_1 = \\gamma_2\$"))
    grid(true)
    legend()
    tight_layout()

    savefig("RaRbGamma2RestrictedPlannerLog.pdf")

end
    
function RaRbConstraint(SOp,Dp)

    # both values set to 0 
    allPlanner = computeAllPlanner(SOp,Dp,RaImpact = 0.0, RbImpact = 0.0)
    allPlannerNoResources = computeAllPlannerNoResources(SOp,Dp,RaImpact = 0.0, RbImpact = 0.0)
    allNash = computeAllNash(SOp,Dp,RaImpact = 0.0, RbImpact = 0.0)
    # allStackelberg = computeAllStackelberg(SOp,Dp)
    keys = [:G,:G1,:G2, :DeltaTemp,:WelFareTot,:WelFare1,:WelFare2]
    TabPlanner = TableCol("Planner", keys, [allPlanner[[11,12,13,16,17,18,19]]...])
    TabPlannerNoResources = TableCol("RestrictedPlanner", keys, [allPlannerNoResources[[11,12,13,16,17,18,19]]...])
    TabNash = TableCol("Nash", keys, [allNash[[11,12,13,16,17,18,19]]...])
    # TabStackelberg = TableCol("Stackelberg", keys, [allStackelberg...])
    # TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash,TabStackelberg)
    TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash)
    to_tex(TableTotal) |> print


    # Ra = 0 
    allPlanner = computeAllPlanner(SOp,Dp,RaImpact = 0.0, RbImpact = 1.0)
    allPlannerNoResources = computeAllPlannerNoResources(SOp,Dp,RaImpact = 0.0, RbImpact = 1.0)
    allNash = computeAllNash(SOp,Dp,RaImpact = 0.0, RbImpact = 1.0)
    # allStackelberg = computeAllStackelberg(SOp,Dp)
    keys = [:G,:G1,:G2, :DeltaTemp,:WelFareTot,:WelFare1,:WelFare2]
    TabPlanner = TableCol("Planner", keys, [allPlanner[[11,12,13,16,17,18,19]]...])
    TabPlannerNoResources = TableCol("RestrictedPlanner", keys, [allPlannerNoResources[[11,12,13,16,17,18,19]]...])
    TabNash = TableCol("Nash", keys, [allNash[[11,12,13,16,17,18,19]]...])
    # TabStackelberg = TableCol("Stackelberg", keys, [allStackelberg...])
    # TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash,TabStackelberg)
    TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash)
    to_tex(TableTotal) |> print

    # Rb = 0 
    allPlanner = computeAllPlanner(SOp,Dp,RaImpact = 1.0, RbImpact = 0.0)
    allPlannerNoResources = computeAllPlannerNoResources(SOp,Dp,RaImpact = 1.0, RbImpact = 0.0)
    allNash = computeAllNash(SOp,Dp,RaImpact = 1.0, RbImpact = 0.0)
    # allStackelberg = computeAllStackelberg(SOp,Dp)
    keys = [:G,:G1,:G2, :DeltaTemp,:WelFareTot,:WelFare1,:WelFare2]
    TabPlanner = TableCol("Planner", keys, [allPlanner[[11,12,13,16,17,18,19]]...])
    TabPlannerNoResources = TableCol("RestrictedPlanner", keys, [allPlannerNoResources[[11,12,13,16,17,18,19]]...])
    TabNash = TableCol("Nash", keys, [allNash[[11,12,13,16,17,18,19]]...])
    # TabStackelberg = TableCol("Stackelberg", keys, [allStackelberg...])
    # TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash,TabStackelberg)
    TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash)
    to_tex(TableTotal) |> print

    # all free
    allPlanner = computeAllPlanner(SOp,Dp,RaImpact = 1.0, RbImpact = 1.0)
    allPlannerNoResources = computeAllPlannerNoResources(SOp,Dp,RaImpact = 1.0, RbImpact = 1.0)
    allNash = computeAllNash(SOp,Dp,RaImpact = 1.0, RbImpact = 1.0)
    # allStackelberg = computeAllStackelberg(SOp,Dp)
    keys = [:G,:G1,:G2, :DeltaTemp,:WelFareTot,:WelFare1,:WelFare2]
    TabPlanner = TableCol("Planner", keys, [allPlanner[[11,12,13,16,17,18,19]]...])
    TabPlannerNoResources = TableCol("RestrictedPlanner", keys, [allPlannerNoResources[[11,12,13,16,17,18,19]]...])
    TabNash = TableCol("Nash", keys, [allNash[[11,12,13,16,17,18,19]]...])
    # TabStackelberg = TableCol("Stackelberg", keys, [allStackelberg...])
    # TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash,TabStackelberg)
    TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash)
    to_tex(TableTotal) |> print

    return nothing
end

function Plannersρ1ρ2(SOp,Dp; Nt = 10)

    times = LinRange(0,Dp.TSim,Nt)


    C1Planner,C2Planner,B1Planner,B2Planner,K1Planner,K2Planner,RaPlanner,RbPlanner = [zeros(Nt) for i in 1:8]
    C1PlannerNoResources,C2PlannerNoResources,B1PlannerNoResources,B2PlannerNoResources,K1PlannerNoResources,K2PlannerNoResources,RaPlannerNoResources,RbPlannerNoResources = [zeros(Nt) for i in 1:8]
    
    @showprogress for (i,t) in enumerate(times)
        Planner_t = optimPlannerDifferntRho(SOp,t)
        C1Planner[i],C2Planner[i],B1Planner[i],B2Planner[i],K1Planner[i],K2Planner[i],RaPlanner[i],RbPlanner[i] = [Planner_t[name] for name in varNames]
        PlannerNoResources_t = optimPlannerNoResourcesDifferentRho(SOp,t)
        C1PlannerNoResources[i],C2PlannerNoResources[i],B1PlannerNoResources[i],B2PlannerNoResources[i],K1PlannerNoResources[i],K2PlannerNoResources[i],RaPlannerNoResources[i],RbPlannerNoResources[i] = [PlannerNoResources_t[name] for name in varNames]
    end

    return times,C1Planner,C2Planner,B1Planner,B2Planner,K1Planner,K2Planner,RaPlanner,RbPlanner,C1PlannerNoResources,C2PlannerNoResources,B1PlannerNoResources,B2PlannerNoResources,K1PlannerNoResources,K2PlannerNoResources,RaPlannerNoResources,RbPlannerNoResources
end

function plotPlannersρ1ρ2(SOp,Dp; Nt = 10)

    times,C1Planner,C2Planner,B1Planner,B2Planner,K1Planner,K2Planner,RaPlanner,RbPlanner,C1PlannerNoResources,C2PlannerNoResources,B1PlannerNoResources,B2PlannerNoResources,K1PlannerNoResources,K2PlannerNoResources,RaPlannerNoResources,RbPlannerNoResources = Plannersρ1ρ2(SOp,Dp,Nt=Nt)
    timesFine = LinRange(times[1],times[end],1000)

    # deleteat!(times,8)
    # deleteat!(C1Planner,8); deleteat!(C1PlannerNoResources,8)
    # deleteat!(C2Planner,8); deleteat!(C2PlannerNoResources,8)
    # deleteat!(B1Planner,8); deleteat!(B1PlannerNoResources,8)
    # deleteat!(B2Planner,8); deleteat!(B2PlannerNoResources,8)
    # deleteat!(K1Planner,8); deleteat!(K1PlannerNoResources,8)
    # deleteat!(K2Planner,8); deleteat!(K2PlannerNoResources,8)
    # deleteat!(RaPlanner,8); deleteat!(RaPlannerNoResources,8)
    # deleteat!(RbPlanner,8); deleteat!(RbPlannerNoResources,8)
    

    # C1
    C1Planner = cubic_spline_interpolation(times,C1Planner)
    C1PlannerNoResources = cubic_spline_interpolation(times,C1PlannerNoResources)
    figure()
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    plot(timesFine, C1Planner.(timesFine), label = "Global Planner - C1",color="blue")
    plot(timesFine, C1PlannerNoResources.(timesFine), label = "Restricted Planner - C1",color="gold")   
    xlabel("Time (years)")
    ylabel("C1")
    grid(true)
    legend()
    tight_layout()
    savefig("C1PlannerGlobalRestrictedDifferentRho.pdf")
   

    # C2
    C2Planner = cubic_spline_interpolation(times,C2Planner)
    C2PlannerNoResources = cubic_spline_interpolation(times,C2PlannerNoResources)
    figure()
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    plot(timesFine, C2Planner.(timesFine), label = "Global Planner - C2",color="blue")
    plot(timesFine, C2PlannerNoResources.(timesFine), label = "Restricted Planner - C2",color="gold")   
    xlabel("Time (years)")
    ylabel("C2")
    grid(true)
    legend()
    tight_layout()
    savefig("C2PlannerGlobalRestrictedDifferentRho.pdf")
   

    # K1
    K1Planner = cubic_spline_interpolation(times,K1Planner)
    K1PlannerNoResources = cubic_spline_interpolation(times,K1PlannerNoResources)
    figure()
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    plot(timesFine, K1Planner.(timesFine), label = "Global Planner - I1",color="blue")
    plot(timesFine, K1PlannerNoResources.(timesFine), label = "Restricted Planner - I1",color="gold")   
    xlabel("Time (years)")
    ylabel("I1")
    grid(true)
    legend()
    tight_layout()
    savefig("I1PlannerGlobalRestrictedDifferentRho.pdf")
   

    # K2
    K2Planner = cubic_spline_interpolation(times,K2Planner)
    K2PlannerNoResources = cubic_spline_interpolation(times,K2PlannerNoResources)
    figure()
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    plot(timesFine, K2Planner.(timesFine), label = "Global Planner - I2",color="blue")
    plot(timesFine, K2PlannerNoResources.(timesFine), label = "Restricted Planner - I2",color="gold")   
    xlabel("Time (years)")
    ylabel("I2")
    grid(true)
    legend()
    tight_layout()
    savefig("I2PlannerGlobalRestrictedDifferentRho.pdf")
   

    # B1
    B1Planner = cubic_spline_interpolation(times,B1Planner)
    B1PlannerNoResources = cubic_spline_interpolation(times,B1PlannerNoResources)
    figure()
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    # plot(timesFine, B1Planner.(timesFine), label = "Global Planner - B1",color="blue")
    plot(timesFine, B1PlannerNoResources.(timesFine), label = "Global Planner - B1",color="blue")
    plot(timesFine, B1PlannerNoResources.(timesFine), label = "Restricted Planner - B1",color="gold")   
    xlabel("Time (years)")
    ylabel("B1")
    grid(true)
    legend()
    tight_layout()
    savefig("B1PlannerGlobalRestrictedDifferentRho.pdf")
   

    # B2
    B2Planner = cubic_spline_interpolation(times,B2Planner)
    B2PlannerNoResources = cubic_spline_interpolation(times,B2PlannerNoResources)
    figure()
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    plot(timesFine, B2Planner.(timesFine), label = "Global Planner - B2",color="blue")
    plot(timesFine, B2PlannerNoResources.(timesFine), label = "Restricted Planner - B2",color="gold")   
    xlabel("Time (years)")
    ylabel("B2")
    grid(true)
    legend()
    tight_layout()
    savefig("B2PlannerGlobalRestrictedDifferentRho.pdf")
   

    # Ra
    RaPlanner = cubic_spline_interpolation(times,RaPlanner)
    RaPlannerNoResources = cubic_spline_interpolation(times,RaPlannerNoResources)
    figure()
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    plot(timesFine, RaPlanner.(timesFine), label = "Global Planner - Ra",color="blue")
    plot(timesFine, RaPlannerNoResources.(timesFine), label = "Restricted Planner - Ra",color="gold")   
    xlabel("Time (years)")
    ylabel("Ra")
    grid(true)
    legend()
    tight_layout()
    savefig("RaPlannerGlobalRestrictedDifferentRho.pdf")
   

    # Rb
    RbPlanner = cubic_spline_interpolation(times,RbPlanner)
    RbPlannerNoResources = cubic_spline_interpolation(times,RbPlannerNoResources)
    figure()
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    plot(timesFine, RbPlanner.(timesFine), label = "Global Planner - Rb",color="blue")
    plot(timesFine, RbPlannerNoResources.(timesFine), label = "Restricted Planner - Rb",color="gold")   
    xlabel("Time (years)")
    ylabel("Rb")
    grid(true)
    legend()
    tight_layout()
    savefig("RbPlannerGlobalRestrictedDifferentRho.pdf")
   
    
    close("all")

end


function optimPlannerDifferntRho(p,t; RaImpact = 1.0, RbImpact = 1.0, quiet = true, showObj = false)
    
    plannerProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerProblem, "tol", GLOBAL_TOL)
    set_optimizer_attribute(plannerProblem, "acceptable_tol",  GLOBAL_TOL)
    set_optimizer_attribute(plannerProblem, "max_iter", GLOBAL_MAX_IT)

    JuMP.@variables(plannerProblem,begin
        C1 ≥ 0.0
        C2 ≥ 0.0
        B1 ≥ 0.0
        B2 ≥ 0.0
        K1 ≥ 0.0
        K2 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    # subexpressions
    @NLexpressions(plannerProblem, begin
        gRb, (p.g∞ * RbImpact*Rb + p.g0) / (RbImpact*Rb + 1)
        hRa, (p.h∞ * RaImpact*Ra + p.h0) / (RaImpact*Ra + 1)
        fK1, K1^p.α
        fK2, K2^p.α
        D1, B1^p.θ1
        D2, gRb * B2^p.θ2
        K, K1 + K2
        D, D1 + D2
        GID, p.ηK * K - p.ηB * D
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
    C1 + C2 + B1 + B2 + K1 + K2 + Ra + Rb ≤ p.A̅ * (fK1 + hRa * fK2))

    # objective
    @NLobjective(plannerProblem, Max, exp(-t*p.ρ1)*(u1C1-p.γ1*p.Φ1*GID) + exp(-t*p.ρ2)*(u2C2-p.γ2*p.Φ2*GID))

    # verbose
    if quiet == false
        unset_silent(plannerProblem)
    else
        set_silent(plannerProblem)
    end

    JuMP.optimize!(plannerProblem)
    
    solValues = value.([C1,C2,B1,B2,K1,K2,Ra,Rb])
    solValues[solValues .< 1e-6] .= 0
    solDict = Dict(zip(varNames,solValues))

    if showObj println("Objective: ",computeWelFarePlanner(p,solDict)) end

    if termination_status(plannerProblem) ∉ [LOCALLY_SOLVED, OPTIMAL]
        printstyled("plannerProblem, termination status -> " * string(termination_status(plannerProblem)) * "\n",color=:red)
        return Dict(zip(varNames,[NaN for i in 1:8]))
    end

    return solDict
end

function optimPlannerNoResourcesDifferentRho(p,t; RaImpact = 1.0, RbImpact = 1.0, quiet = true, start = nothing, showObj = false)
    
    plannerNRProblem = Model(Ipopt.Optimizer)
    set_optimizer_attribute(plannerNRProblem, "tol", GLOBAL_TOL)
    set_optimizer_attribute(plannerNRProblem, "acceptable_tol",  GLOBAL_TOL)
    set_optimizer_attribute(plannerNRProblem, "max_iter", GLOBAL_MAX_IT)

    JuMP.@variables(plannerNRProblem,begin
        C1 ≥ 0.0
        C2 ≥ 0.0
        B1 ≥ 0.0
        B2 ≥ 0.0
        K1 ≥ 0.0
        K2 ≥ 0.0
        Ra ≥ 0.0
        Rb ≥ 0.0
    end)

    if !isnothing(start)
        set_start_value(C1,start["C1"])
        set_start_value(C2,start["C2"])
        set_start_value(B1,start["B1"])
        set_start_value(B2,start["B2"])
        set_start_value(K1,start["K1"])
        set_start_value(K2,start["K2"])
        set_start_value(Ra,start["Ra"])
        set_start_value(Rb,start["Rb"])
    end

    # subexpressions
    @NLexpressions(plannerNRProblem, begin
    gRb, (p.g∞ * RbImpact*Rb + p.g0) / (RbImpact*Rb + 1)
    hRa, (p.h∞ * RaImpact*Ra + p.h0) / (RaImpact*Ra + 1)
    fK1, K1^p.α
    fK2, K2^p.α
    D1, B1^p.θ1
    D2, gRb * B2^p.θ2
    K, K1 + K2
    D, D1 + D2
    GID, p.ηK * K - p.ηB * D
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
    C1 + B1 + K1 + Ra + Rb ≤ p.A̅ * fK1)

    # budget constraint 2
    @NLconstraint(plannerNRProblem,BudgetConstraint2,
    C2 + B2 + K2 ≤ p.A̅ * hRa * fK2)
    
    # objective
    @NLobjective(plannerNRProblem, Max, exp(-t*p.ρ1)*(u1C1-p.γ1*p.Φ1*GID) + exp(-t*p.ρ2)*(u2C2-p.γ2*p.Φ2*GID))

    # verbose
    if quiet == false
        unset_silent(plannerNRProblem)
    else
        set_silent(plannerNRProblem)
    end

    JuMP.optimize!(plannerNRProblem)
        
    solValues = value.([C1,C2,B1,B2,K1,K2,Ra,Rb])
    solValues[solValues .< 1e-6] .= 0
    solDict = Dict(zip(varNames,solValues))

    if showObj println("Objective: ",computeWelFarePlannerNoResources(p,solDict)) end

    if termination_status(plannerNRProblem) ∉ [LOCALLY_SOLVED, OPTIMAL]
        printstyled("plannerNRProblem, termination status -> " * string(termination_status(plannerNRProblem)) * "\n",color=:red)
        return Dict(zip(varNames,[NaN for i in 1:8]))
    end
    
    return solDict
end