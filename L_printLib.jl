
function makeTableResultsExplicit(SOp,Dp)

    allPlanner = computeAllPlannerExplicit(SOp,Dp)
    allPlannerNoResources = computeAllPlannerNoResourcesExplicit(SOp,Dp)
    allNash = computeAllNashExplicit(SOp,Dp)
    allStackelberg = computeAllStackelbergExplicit(SOp,Dp)

    keys = [[Symbol(name) for name in varNames]...,:gRb,:hRa,:G,:G1,:G2, :DeltaP,:DeltaT,:DeltaTemp,:WelFareTot,:WelFare1,:WelFare2,:Y1,:Y2]

    TabPlanner = TableCol("Planner", keys, [allPlanner...])
    TabPlannerNoResources = TableCol("RestrictedPlanner", keys, [allPlannerNoResources...])
    TabNash = TableCol("Nash", keys, [allNash...])
    TabStackelberg = TableCol("Stackelberg", keys, [allStackelberg...])

    TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash,TabStackelberg)
    to_tex(TableTotal) |> print

    return nothing
end

function makeTableResults(SOp,Dp)

    allPlanner = computeAllPlanner(SOp,Dp)
    allPlannerNoResources = computeAllPlannerNoResources(SOp,Dp)
    allNash = computeAllNash(SOp,Dp)
    allStackelberg = computeAllStackelberg(SOp,Dp)

    keys = [[Symbol(name) for name in varNames]...,:G, :DeltaP,:DeltaT,:DeltaTemp,:WelFareTot,:WelFare1,:WelFare2,:Y1,:Y2]

    TabPlanner = TableCol("Planner", keys, [allPlanner...])
    TabPlannerNoResources = TableCol("RestrictedPlanner", keys, [allPlannerNoResources...])
    TabNash = TableCol("Nash", keys, [allNash...])
    TabStackelberg = TableCol("Stackelberg", keys, [allStackelberg...])

    TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash,TabStackelberg)
    to_tex(TableTotal) |> print

    return nothing
end

function computeTemperatureAllExplicit(SOp,Dp)

    solPlanner = optimPlannerExplicit(SOp,quiet=true)
    GComputed = computeG(SOp,solPlanner)
    solODEPlanner = solveODE(SOp,Dp,GComputed)
    TempPlanner = [ComputeTemperature(Dp,p,t) for (p,t) in solODEPlanner.u]

    solPlannerNoResources = optimPlannerNoResourcesExplicit(SOp,quiet=true)
    GComputed = computeG(SOp,solPlannerNoResources)
    solODEPlannerNoResources = solveODE(SOp,Dp,GComputed)
    TempPlannerNoResources = [ComputeTemperature(Dp,p,t) for (p,t) in solODEPlannerNoResources.u]

    solNash = optimNashExplicit(SOp,quiet=true)
    GComputed = computeG(SOp,solNash)
    solODENash = solveODE(SOp,Dp,GComputed)
    TempNash = [ComputeTemperature(Dp,p,t) for (p,t) in solODENash.u]

    solStackelberg = optimStackelbergExplicit(SOp,quiet=true)
    GComputed = computeG(SOp,solStackelberg)
    solODEStackelberg = solveODE(SOp,Dp,GComputed)
    TempStackelberg = [ComputeTemperature(Dp,p,t) for (p,t) in solODEStackelberg.u]

    return TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg
end


function computeTemperatureAll(SOp,Dp)

    solPlanner = optimPlanner(SOp,quiet=true)
    GComputed = computeG(SOp,solPlanner)
    solODEPlanner = solveODE(SOp,Dp,GComputed)
    TempPlanner = [ComputeTemperature(Dp,p,t) for (p,t) in solODEPlanner.u]

    solPlannerNoResources = optimPlannerNoResources(SOp,quiet=true)
    GComputed = computeG(SOp,solPlannerNoResources)
    solODEPlannerNoResources = solveODE(SOp,Dp,GComputed)
    TempPlannerNoResources = [ComputeTemperature(Dp,p,t) for (p,t) in solODEPlannerNoResources.u]

    solNash = optimNash(SOp,quiet=true)
    GComputed = computeG(SOp,solNash)
    solODENash = solveODE(SOp,Dp,GComputed)
    TempNash = [ComputeTemperature(Dp,p,t) for (p,t) in solODENash.u]

    solStackelberg = optimStackelberg(SOp,quiet=true)
    GComputed = computeG(SOp,solStackelberg)
    solODEStackelberg = solveODE(SOp,Dp,GComputed)
    TempStackelberg = [ComputeTemperature(Dp,p,t) for (p,t) in solODEStackelberg.u]

    return TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg
end

function computeTemperatureAllRobust(SOp,Dp; N = 100, quiet = true)

    SOpAux = deepcopy(SOp)
    XExp = Exponential(1)
    TempPlannerRobust = zeros(N,length(Dp.TSave))
    TempPlannerNoResourcesRobust = zeros(N,length(Dp.TSave))
    TempNashRobust = zeros(N,length(Dp.TSave))
    TempStackelbergRobust = zeros(N,length(Dp.TSave))

    @showprogress for i in 1:N
        γ1sample = SOp.γ1 * rand(XExp)
        γ2sample = SOp.γ2 * rand(XExp)
        
        SOpAux.γ̂1 = γ1sample; SOpAux.γ1 = γ1sample; SOpAux.γ̂2 = γ2sample; SOpAux.γ2 = γ2sample
        solPlannerRobust, γ1Robust, γ2Robust = optimPlannerRobust(SOpAux,Dp; quiet = quiet)
        SOpAux.γ̂1 = γ1Robust; SOpAux.γ1 = γ1Robust; SOpAux.γ̂2 = γ2Robust; SOpAux.γ2 = γ2Robust
        GComputed = computeG(SOpAux,solPlannerRobust)
        solODEPlannerRobust = solveODE(SOpAux,Dp,GComputed)
        TempPlannerRobust[i,:] = [ComputeTemperature(Dp,p,t) for (p,t) in solODEPlannerRobust.u]

        SOpAux.γ̂1 = γ1sample; SOpAux.γ1 = γ1sample; SOpAux.γ̂2 = γ2sample; SOpAux.γ2 = γ2sample
        solPlannerNoResourcesRobust, γ1Robust, γ2Robust = optimPlannerNoResourcesRobust(SOpAux,Dp; quiet = quiet)
        SOpAux.γ̂1 = γ1Robust; SOpAux.γ1 = γ1Robust; SOpAux.γ̂2 = γ2Robust; SOpAux.γ2 = γ2Robust
        GComputed = computeG(SOpAux,solPlannerNoResourcesRobust)
        solODEPlannerNoResourcesRobust = solveODE(SOpAux,Dp,GComputed)
        TempPlannerNoResourcesRobust[i,:] = [ComputeTemperature(Dp,p,t) for (p,t) in solODEPlannerNoResourcesRobust.u]
    
        SOpAux.γ̂1 = γ1sample; SOpAux.γ1 = γ1sample; SOpAux.γ̂2 = γ2sample; SOpAux.γ2 = γ2sample
        solNashRobust, γ1Robust, γ2Robust = optimNashRobust(SOpAux,Dp; quiet = quiet)
        SOpAux.γ̂1 = γ1Robust; SOpAux.γ1 = γ1Robust; SOpAux.γ̂2 = γ2Robust; SOpAux.γ2 = γ2Robust
        GComputed = computeG(SOpAux,solNashRobust)
        solODENashRobust = solveODE(SOpAux,Dp,GComputed)
        TempNashRobust[i,:] = [ComputeTemperature(Dp,p,t) for (p,t) in solODENashRobust.u]
    
        SOpAux.γ̂1 = γ1sample; SOpAux.γ1 = γ1sample; SOpAux.γ̂2 = γ2sample; SOpAux.γ2 = γ2sample
        solStackelbergRobust, γ1Robust, γ2Robust = optimStackelbergRobust(SOpAux,Dp; quiet = quiet)
        SOpAux.γ̂1 = γ1Robust; SOpAux.γ1 = γ1Robust; SOpAux.γ̂2 = γ2Robust; SOpAux.γ2 = γ2Robust
        GComputed = computeG(SOpAux,solStackelbergRobust)
        solODEStackelbergRobust = solveODE(SOpAux,Dp,GComputed)
        TempStackelbergRobust[i,:] = [ComputeTemperature(Dp,p,t) for (p,t) in solODEStackelbergRobust.u]
    end


    return TempPlannerRobust, TempPlannerNoResourcesRobust, TempNashRobust, TempStackelbergRobust
end


function computeEmissionAllExplicit(SOp,Dp)

    solPlanner = optimPlannerExplicit(SOp,quiet=true)
    GComputed = computeG(SOp,solPlanner)
    G1Planner = computeG1(SOp,solPlanner)
    G2Planner = computeG2(SOp,solPlanner)
    Dp1Planner,Dp2Planner = splitInitialData(Dp,G1Planner,G2Planner)
    solODEPlanner1 = solveODE(SOp,Dp1Planner,G1Planner)
    solODEPlanner2 = solveODE(SOp,Dp2Planner,G2Planner)
    pPlanner1, tPlanner1 = solODEPlanner1[1,:], solODEPlanner1[2,:]
    pPlanner2, tPlanner2 = solODEPlanner2[1,:], solODEPlanner2[2,:]
    pPlanner = pPlanner1 + pPlanner2
    tPlanner = tPlanner1 + tPlanner2

    solPlannerNoResources = optimPlannerNoResourcesExplicit(SOp,quiet=true)
    GComputed = computeG(SOp,solPlannerNoResources)
    G1PlannerNoResources = computeG1(SOp,solPlannerNoResources)
    G2PlannerNoResources = computeG2(SOp,solPlannerNoResources)
    Dp1PlannerNoResources,Dp2PlannerNoResources = splitInitialData(Dp,G1PlannerNoResources,G2PlannerNoResources)
    solODEPlannerNoResources1 = solveODE(SOp,Dp1PlannerNoResources,G1PlannerNoResources)
    solODEPlannerNoResources2 = solveODE(SOp,Dp2PlannerNoResources,G2PlannerNoResources)
    pPlannerNoResources1, tPlannerNoResources1 = solODEPlannerNoResources1[1,:], solODEPlannerNoResources1[2,:]
    pPlannerNoResources2, tPlannerNoResources2 = solODEPlannerNoResources2[1,:], solODEPlannerNoResources2[2,:]
    pPlannerNoResources = pPlannerNoResources1 + pPlannerNoResources2
    tPlannerNoResources = tPlannerNoResources1 + tPlannerNoResources2

    solNash = optimNashExplicit(SOp,quiet=true)
    GComputed = computeG(SOp,solNash)
    G1Nash = computeG1(SOp,solNash)
    G2Nash = computeG2(SOp,solNash)
    Dp1Nash,Dp2Nash = splitInitialData(Dp,G1Nash,G2Nash)
    solODENash1 = solveODE(SOp,Dp1Nash,G1Nash)
    solODENash2 = solveODE(SOp,Dp2Nash,G2Nash)
    pNash1, tNash1 = solODENash1[1,:], solODENash1[2,:]
    pNash2, tNash2 = solODENash2[1,:], solODENash2[2,:]
    pNash = pNash1 + pNash2
    tNash = tNash1 + tNash2

    solStackelberg = optimStackelbergExplicit(SOp,quiet=true)
    GComputed = computeG(SOp,solStackelberg)
    G1Stackelberg = computeG1(SOp,solStackelberg)
    G2Stackelberg = computeG2(SOp,solStackelberg)
    Dp1Stackelberg,Dp2Stackelberg = splitInitialData(Dp,G1Stackelberg,G2Stackelberg)
    solODEStackelberg1 = solveODE(SOp,Dp1Stackelberg,G1Stackelberg)
    solODEStackelberg2 = solveODE(SOp,Dp2Stackelberg,G2Stackelberg)
    pStackelberg1, tStackelberg1 = solODEStackelberg1[1,:], solODEStackelberg1[2,:]
    pStackelberg2, tStackelberg2 = solODEStackelberg2[1,:], solODEStackelberg2[2,:]
    pStackelberg = pStackelberg1 + pStackelberg2
    tStackelberg = tStackelberg1 + tStackelberg2

    return pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2

end



function computeEmissionAll(SOp,Dp)

    solPlanner = optimPlanner(SOp,quiet=true)
    GComputed = computeG(SOp,solPlanner)
    G1Planner = computeG1(SOp,solPlanner)
    G2Planner = computeG2(SOp,solPlanner)
    Dp1Planner,Dp2Planner = splitInitialData(Dp,G1Planner,G2Planner)
    solODEPlanner1 = solveODE(SOp,Dp1Planner,G1Planner)
    solODEPlanner2 = solveODE(SOp,Dp2Planner,G2Planner)
    pPlanner1, tPlanner1 = solODEPlanner1[1,:], solODEPlanner1[2,:]
    pPlanner2, tPlanner2 = solODEPlanner2[1,:], solODEPlanner2[2,:]
    pPlanner = pPlanner1 + pPlanner2
    tPlanner = tPlanner1 + tPlanner2

    solPlannerNoResources = optimPlannerNoResources(SOp,quiet=true)
    GComputed = computeG(SOp,solPlannerNoResources)
    G1PlannerNoResources = computeG1(SOp,solPlannerNoResources)
    G2PlannerNoResources = computeG2(SOp,solPlannerNoResources)
    Dp1PlannerNoResources,Dp2PlannerNoResources = splitInitialData(Dp,G1PlannerNoResources,G2PlannerNoResources)
    solODEPlannerNoResources1 = solveODE(SOp,Dp1PlannerNoResources,G1PlannerNoResources)
    solODEPlannerNoResources2 = solveODE(SOp,Dp2PlannerNoResources,G2PlannerNoResources)
    pPlannerNoResources1, tPlannerNoResources1 = solODEPlannerNoResources1[1,:], solODEPlannerNoResources1[2,:]
    pPlannerNoResources2, tPlannerNoResources2 = solODEPlannerNoResources2[1,:], solODEPlannerNoResources2[2,:]
    pPlannerNoResources = pPlannerNoResources1 + pPlannerNoResources2
    tPlannerNoResources = tPlannerNoResources1 + tPlannerNoResources2

    solNash = optimNash(SOp,quiet=true)
    GComputed = computeG(SOp,solNash)
    G1Nash = computeG1(SOp,solNash)
    G2Nash = computeG2(SOp,solNash)
    Dp1Nash,Dp2Nash = splitInitialData(Dp,G1Nash,G2Nash)
    solODENash1 = solveODE(SOp,Dp1Nash,G1Nash)
    solODENash2 = solveODE(SOp,Dp2Nash,G2Nash)
    pNash1, tNash1 = solODENash1[1,:], solODENash1[2,:]
    pNash2, tNash2 = solODENash2[1,:], solODENash2[2,:]
    pNash = pNash1 + pNash2
    tNash = tNash1 + tNash2

    solStackelberg = optimStackelberg(SOp,quiet=true)
    GComputed = computeG(SOp,solStackelberg)
    G1Stackelberg = computeG1(SOp,solStackelberg)
    G2Stackelberg = computeG2(SOp,solStackelberg)
    Dp1Stackelberg,Dp2Stackelberg = splitInitialData(Dp,G1Stackelberg,G2Stackelberg)
    solODEStackelberg1 = solveODE(SOp,Dp1Stackelberg,G1Stackelberg)
    solODEStackelberg2 = solveODE(SOp,Dp2Stackelberg,G2Stackelberg)
    pStackelberg1, tStackelberg1 = solODEStackelberg1[1,:], solODEStackelberg1[2,:]
    pStackelberg2, tStackelberg2 = solODEStackelberg2[1,:], solODEStackelberg2[2,:]
    pStackelberg = pStackelberg1 + pStackelberg2
    tStackelberg = tStackelberg1 + tStackelberg2

    return pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2

end



function plotTemperature(SOp,Dp;explicit=true)
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    if explicit & (SOp.α == 1)
        TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg = computeTemperatureAllExplicit(SOp,Dp)        
    else
        TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg = computeTemperatureAll(SOp,Dp)
    end

    plot(Dp.TSave, TempPlanner, label = "Global planner",color="blue")
    plot(Dp.TSave, TempPlannerNoResources, label = "Restricted planner",color="gold")
    plot(Dp.TSave, TempNash, label = "Nash",color="green")
    plot(Dp.TSave, TempStackelberg, label = "Stackelberg",color="red")
    xlabel("Time (years)")
    ylabel("Global temperature (°C)")
    grid()
    legend()
    tight_layout()

    savefig("Temperature.pdf")
end

function plotTemperatureRobustPoint(SOp,Dp; quiet = true)
    alphaRange = [10^i for i in 1:6]
    # alphaRange = [10.0^i for i in [1,2,3,4]]

    N = length(alphaRange)
    SOpAux = deepcopy(SOp)
   
    TempPlannerRobust = zeros(N,length(Dp.TSave))
    TempPlannerNoResourcesRobust = zeros(N,length(Dp.TSave))
    TempNashRobust = zeros(N,length(Dp.TSave))
    TempStackelbergRobust = zeros(N,length(Dp.TSave))

    for i in 1:N
        
        SOpAux.αR = alphaRange[i]

        solPlannerRobust, γ1Robust, γ2Robust = optimPlannerRobust(SOpAux,Dp; quiet = quiet)
        GComputed = computeG(SOpAux,solPlannerRobust)
        solODEPlannerRobust = solveODE(SOpAux,Dp,GComputed)
        TempPlannerRobust[i,:] = [ComputeTemperature(Dp,p,t) for (p,t) in solODEPlannerRobust.u]

        solPlannerNoResourcesRobust, γ1Robust, γ2Robust = optimPlannerNoResourcesRobust(SOpAux,Dp; quiet = quiet)
        GComputed = computeG(SOpAux,solPlannerNoResourcesRobust)
        solODEPlannerNoResourcesRobust = solveODE(SOpAux,Dp,GComputed)
        TempPlannerNoResourcesRobust[i,:] = [ComputeTemperature(Dp,p,t) for (p,t) in solODEPlannerNoResourcesRobust.u]
    
        solNashRobust, γ1Robust, γ2Robust = optimNashRobust(SOpAux,Dp; quiet = quiet)
        GComputed = computeG(SOpAux,solNashRobust)
        solODENashRobust = solveODE(SOpAux,Dp,GComputed)
        TempNashRobust[i,:] = [ComputeTemperature(Dp,p,t) for (p,t) in solODENashRobust.u]
    
        solStackelbergRobust, γ1Robust, γ2Robust = optimStackelbergRobust(SOpAux,Dp; quiet = quiet)
        GComputed = computeG(SOpAux,solStackelbergRobust)
        solODEStackelbergRobust = solveODE(SOpAux,Dp,GComputed)
        TempStackelbergRobust[i,:] = [ComputeTemperature(Dp,p,t) for (p,t) in solODEStackelbergRobust.u]
    end

    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    for i in 1:N
        if i == N
            label = "Global planner, " * "α = ∞"
        else            
            label = "Global planner, " * "α = " * string(round(alphaRange[i],digits=2)) 
        end
        plot(Dp.TSave, TempPlannerRobust[i,:], label = label,color="blue", ls = ALL_LINE_STYLES[i])
    end
    xlabel("Time (years)")
    ylabel("Global temperature (°C)")
    grid(true)
    legend()
    tight_layout()
    savefig("TemperatureRobustPointPlanner.pdf")
    fig.clear()

    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    for i in 1:N
        if i == N
            label = "Restricted planner, " * "α = ∞"
        else            
            label = "Restricted planner, " * "α = " * string(round(alphaRange[i],digits=2)) 
        end
        plot(Dp.TSave, TempPlannerNoResourcesRobust[i,:], label = label,color="gold", ls = ALL_LINE_STYLES[i])
    end
    xlabel("Time (years)")
    ylabel("Global temperature (°C)")
    grid(true)
    legend()
    tight_layout()
    savefig("TemperatureRobustPointPlannerNoResources.pdf")
    fig.clear()

    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    for i in 1:N
        if i == N
            label = "Nash, " * "α = ∞"
        else            
            label = "Nash, " * "α = " * string(round(alphaRange[i],digits=2)) 
        end
        plot(Dp.TSave, TempNashRobust[i,:], label = label,color="green", ls = ALL_LINE_STYLES[i])
    end
    xlabel("Time (years)")
    ylabel("Global temperature (°C)")
    grid(true)
    legend()
    tight_layout()
    savefig("TemperatureRobustPointNash.pdf")
    fig.clear()

    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    for i in 1:N
        if i == N
            label = "Stackelberg, " * "α = ∞"
        else            
            label = "Stackelberg, " * "α = " * string(round(alphaRange[i],digits=2)) 
        end
        plot(Dp.TSave, TempStackelbergRobust[i,:], label = label,color="red", ls = ALL_LINE_STYLES[i])
    end
    xlabel("Time (years)")
    ylabel("Global temperature (°C)")
    grid(true)
    legend()
    tight_layout()
    savefig("TemperatureRobustPointStackelberg.pdf")
    fig.clear()


end

function plotTemperatureRobust(SOp,Dp; alpha = 0.05, quiet = true)
    TempPlannerRobust, TempPlannerNoResourcesRobust, TempNashRobust, TempStackelbergRobust = computeTemperatureAllRobust(SOp,Dp, quiet = quiet)
    
    TempPlannerMean = [mean(t) for t in eachcol(TempPlannerRobust)]
    TempPlannerUp = [quantile(t,1-alpha/2) for t in eachcol(TempPlannerRobust)]
    TempPlannerDown = [quantile(t,alpha/2) for t in eachcol(TempPlannerRobust)]
    
    TempPlannerNoResourcesMean = [mean(t) for t in eachcol(TempPlannerNoResourcesRobust)]
    TempPlannerNoResourcesUp = [quantile(t,1-alpha/2) for t in eachcol(TempPlannerNoResourcesRobust)]
    TempPlannerNoResourcesDown = [quantile(t,alpha/2) for t in eachcol(TempPlannerNoResourcesRobust)]

    TempNashMean = [mean(t) for t in eachcol(TempNashRobust)]
    TempNashUp = [quantile(t,1-alpha/2) for t in eachcol(TempNashRobust)]
    TempNashDown = [quantile(t,alpha/2) for t in eachcol(TempNashRobust)]

    TempStackelbergMean = [mean(t) for t in eachcol(TempStackelbergRobust)]
    TempStackelbergUp = [quantile(t,1-alpha/2) for t in eachcol(TempStackelbergRobust)]
    TempStackelbergDown = [quantile(t,alpha/2) for t in eachcol(TempStackelbergRobust)]


    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    fill_between(Dp.TSave, TempPlannerDown, TempPlannerUp, alpha = 0.2,color="blue")
    plot(Dp.TSave, TempPlannerMean, label = "Global planner",color="blue")
    fill_between(Dp.TSave, TempPlannerNoResourcesDown, TempPlannerNoResourcesUp,alpha = 0.2,color="gold")
    plot(Dp.TSave, TempPlannerNoResourcesMean,  label = "Restricted planner",color="gold")
    fill_between(Dp.TSave, TempNashDown, TempNashUp,alpha = 0.2,color="green")
    plot(Dp.TSave, TempNashMean, label = "Nash",color="green")
    fill_between(Dp.TSave, TempStackelbergDown, TempStackelbergUp, alpha = 0.2,color="red")
    plot(Dp.TSave, TempStackelbergMean, label = "Stackelberg",color="red")
    xlabel("Time (years)")
    ylabel("Global temperature (°C)")
    grid()
    legend()
    tight_layout()
    
    savefig("TemperatureRobust.pdf")
end 

function plotTotalEmission(SOp,Dp;explicit = true)
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    if explicit & (SOp.α == 1)
        pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2 = computeEmissionAllExplicit(SOp,Dp)
    else
        pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2 = computeEmissionAll(SOp,Dp)
    end

    # pPlanner1 = pPlanner1 .- pPlanner1[1] .+ pPlanner[1]
    # tPlanner1 = tPlanner1 .- tPlanner1[1] .+ tPlanner[1]
    # pPlanner2 = pPlanner2 .- pPlanner2[1] .+ pPlanner[1]
    # tPlanner2 = tPlanner2 .- tPlanner2[1] .+ tPlanner[1]

    # pPlannerNoResources1 = pPlannerNoResources1 .- pPlannerNoResources1[1] .+ pPlannerNoResources[1]
    # tPlannerNoResources1 = tPlannerNoResources1 .- tPlannerNoResources1[1] .+ tPlannerNoResources[1]
    # pPlannerNoResources2 = pPlannerNoResources2 .- pPlannerNoResources2[1] .+ pPlannerNoResources[1]
    # tPlannerNoResources2 = tPlannerNoResources2 .- tPlannerNoResources2[1] .+ tPlannerNoResources[1]

    # pNash1 = pNash1 .- pNash1[1] .+ pNash[1]
    # tNash1 = tNash1 .- tNash1[1] .+ tNash[1]
    # pNash2 = pNash2 .- pNash2[1] .+ pNash[1]
    # tNash2 = tNash2 .- tNash2[1] .+ tNash[1]

    # pStackelberg1 = pStackelberg1 .- pStackelberg1[1] .+ pStackelberg[1]
    # tStackelberg1 = tStackelberg1 .- tStackelberg1[1] .+ tStackelberg[1]
    # pStackelberg2 = pStackelberg2 .- pStackelberg2[1] .+ pStackelberg[1]
    # tStackelberg2 = tStackelberg2 .- tStackelberg2[1] .+ tStackelberg[1]

    plot(Dp.TSave, pPlanner + tPlanner, label = "Global planner - total emissions",color="blue")
    plot(Dp.TSave, pPlanner1 + tPlanner1, label = "Global planner - country 1 emissions",color="blue", "--")
    plot(Dp.TSave, pPlanner2 + tPlanner2, label = "Global planner - country 2 emissions",color="blue", ":")
    plot(Dp.TSave, pPlannerNoResources + tPlannerNoResources, label = "Restricted planner - total emissions",color="gold")
    plot(Dp.TSave, pPlannerNoResources1 + tPlannerNoResources1, label = "Restricted planner - country 1 emissions",color="gold","--")
    plot(Dp.TSave, pPlannerNoResources2 + tPlannerNoResources2, label = "Restricted planner - country 2 emissions",color="gold",":")
    plot(Dp.TSave, pNash + tNash, label = "Nash - total emissions",color="green")
    plot(Dp.TSave, pNash1 + tNash1, label = "Nash - country 1 emissions",color="green","--")
    plot(Dp.TSave, pNash2 + tNash2, label = "Nash - country 2 emissions",color="green",":")
    plot(Dp.TSave, pStackelberg + tStackelberg, label = "Stackelberg - total emissions",color="red")
    plot(Dp.TSave, pStackelberg1 + tStackelberg1, label = "Stackelberg - country 1 emissions",color="red","--")
    plot(Dp.TSave, pStackelberg2 + tStackelberg2, label = "Stackelberg - country 2 emissions",color="red",":")

    xlabel("Time (years)")
    ylabel("CO2 concentration (ppm)")
    grid()
    legend()
    tight_layout()

    savefig("TotalEmissions.pdf")
end

function plotPTEmission(SOp,Dp;explicit = true)
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    if explicit & (SOp.α == 1)
        pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2 = computeEmissionAllExplicit(SOp,Dp)
    else
        pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2 = computeEmissionAll(SOp,Dp)
    end

    plot(Dp.TSave, pPlanner, label = "Global planner - permanent emission",color="blue","-")
    plot(Dp.TSave, pPlannerNoResources, label = "Restricted planner - permanent emission",color="gold","-")
    plot(Dp.TSave, pNash, label = "Nash - permanent emission",color="green","-")
    plot(Dp.TSave, pStackelberg, label = "Stackelberg - permanent emission",color="red","-")
    legend(loc=7, title="(Left scale)")    
    grid()
    TemporaryPlot = twinx()
    TemporaryPlot.tick_params(axis="y",color="red",labelcolor="red")
    
    plot(Dp.TSave, tPlanner, label = "Global planner - temporary emission",color="blue","--")
    plot(Dp.TSave, tPlannerNoResources, label = "Restricted planner - temporary emission",color="gold","--")
    plot(Dp.TSave, tNash, label = "Nash - temporary emission",color="green","--")
    plot(Dp.TSave, tStackelberg, label = "Stackelberg - temporary emission",color="red","--")

    xlabel("Time (years)")
    ylabel("CO2 concentration (ppm)")

    legend(loc=9,title="(Right scale)") 
    tight_layout()

    savefig("PermTempEmissions.pdf")

end



function plotgh(SOp,Dp;explicit = true)

    if explicit & (SOp.α == 1)
        allPlanner = computeAllPlannerExplicit(SOp,Dp)
        allPlannerNoResources = computeAllPlannerNoResourcesExplicit(SOp,Dp)
        allNash = computeAllNashExplicit(SOp,Dp)
        allStackelberg = computeAllStackelbergExplicit(SOp,Dp)
    else
        allPlanner = computeAllPlanner(SOp,Dp)
        allPlannerNoResources = computeAllPlannerNoResources(SOp,Dp)
        allNash = computeAllNash(SOp,Dp)
        allStackelberg = computeAllStackelberg(SOp,Dp)
    end

    RbPlanner = allPlanner[8]
    RbPlannerNoResources = allPlannerNoResources[8]
    RbNash = allNash[8]
    RbStackelberg = allStackelberg[8]

    RaPlanner = allPlanner[7]
    RaPlannerNoResources = allPlannerNoResources[7]
    RaNash = allNash[7]
    RaStackelberg = allStackelberg[7]

    RabRange = LinRange(1,1.1*max(RbPlanner,RbPlannerNoResources,RbNash,RbStackelberg,RaPlanner,RaPlannerNoResources,RaNash,RaStackelberg),1000)
    semilogx(RabRange,[g(x,SOpG) for x in RabRange],label="g(x)")
    semilogx(RabRange,[h(x,SOpG) for x in RabRange],label="h(x)")
    xlabel("x (log)")
    ylabel(" ")
    grid()
    legend()
    tight_layout()
    savefig("gh.pdf")

end

