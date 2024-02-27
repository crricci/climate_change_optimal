
function makeTableResultsExplicit()

    allPlanner = computeAllPlannerExplicit()
    allPlannerNoResources = computeAllPlannerNoResourcesExplicit()
    allNash = computeAllNashExplicit()
    allStackelberg = computeAllStackelbergExplicit()

    keys = [[Symbol(name) for name in varNames]...,:gRb,:hRa,:G,:G1,:G2, :DeltaP,:DeltaT,:DeltaTemp,:WelFareTot,:WelFare1,:WelFare2,:Y1,:Y2]

    TabPlanner = TableCol("Planner", keys, [allPlanner...])
    TabPlannerNoResources = TableCol("RestrictedPlanner", keys, [allPlannerNoResources...])
    TabNash = TableCol("Nash", keys, [allNash...])
    TabStackelberg = TableCol("Stackelberg", keys, [allStackelberg...])

    TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash,TabStackelberg)
    to_tex(TableTotal) |> print

    return nothing
end

function makeTableResults()

    allPlanner = computeAllPlanner()
    allPlannerNoResources = computeAllPlannerNoResources()
    allNash = computeAllNash()
    allStackelberg = computeAllStackelberg()

    keys = [[Symbol(name) for name in varNames]...,:G, :DeltaP,:DeltaT,:DeltaTemp,:WelFareTot,:WelFare1,:WelFare2,:Y1,:Y2]

    TabPlanner = TableCol("Planner", keys, [allPlanner...])
    TabPlannerNoResources = TableCol("RestrictedPlanner", keys, [allPlannerNoResources...])
    TabNash = TableCol("Nash", keys, [allNash...])
    TabStackelberg = TableCol("Stackelberg", keys, [allStackelberg...])

    TableTotal = hcat(TabPlanner,TabPlannerNoResources,TabNash,TabStackelberg)
    to_tex(TableTotal) |> print

    return nothing
end

function computeTemperatureAllExplicit()

    solPlanner = optimPlannerExplicit(SOpG,quiet=true)
    GComputed = computeG(SOpG,solPlanner)
    solODEPlanner = solveODE(SOpG,DpG,GComputed)
    TempPlanner = [ComputeTemperature(DpG,p,t) for (p,t) in solODEPlanner.u]

    solPlannerNoResources = optimPlannerNoResourcesExplicit(SOpG,quiet=true)
    GComputed = computeG(SOpG,solPlannerNoResources)
    solODEPlannerNoResources = solveODE(SOpG,DpG,GComputed)
    TempPlannerNoResources = [ComputeTemperature(DpG,p,t) for (p,t) in solODEPlannerNoResources.u]

    solNash = optimNashExplicit(SOpG,quiet=true)
    GComputed = computeG(SOpG,solNash)
    solODENash = solveODE(SOpG,DpG,GComputed)
    TempNash = [ComputeTemperature(DpG,p,t) for (p,t) in solODENash.u]

    solStackelberg = optimStackelbergExplicit(SOpG,quiet=true)
    GComputed = computeG(SOpG,solStackelberg)
    solODEStackelberg = solveODE(SOpG,DpG,GComputed)
    TempStackelberg = [ComputeTemperature(DpG,p,t) for (p,t) in solODEStackelberg.u]

    return TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg
end


function computeTemperatureAll()

    solPlanner = optimPlanner(SOpG,quiet=true)
    GComputed = computeG(SOpG,solPlanner)
    solODEPlanner = solveODE(SOpG,DpG,GComputed)
    TempPlanner = [ComputeTemperature(DpG,p,t) for (p,t) in solODEPlanner.u]

    solPlannerNoResources = optimPlannerNoResources(SOpG,quiet=true)
    GComputed = computeG(SOpG,solPlannerNoResources)
    solODEPlannerNoResources = solveODE(SOpG,DpG,GComputed)
    TempPlannerNoResources = [ComputeTemperature(DpG,p,t) for (p,t) in solODEPlannerNoResources.u]

    solNash = optimNash(SOpG,quiet=true)
    GComputed = computeG(SOpG,solNash)
    solODENash = solveODE(SOpG,DpG,GComputed)
    TempNash = [ComputeTemperature(DpG,p,t) for (p,t) in solODENash.u]

    solStackelberg = optimStackelberg(SOpG,quiet=true)
    GComputed = computeG(SOpG,solStackelberg)
    solODEStackelberg = solveODE(SOpG,DpG,GComputed)
    TempStackelberg = [ComputeTemperature(DpG,p,t) for (p,t) in solODEStackelberg.u]

    return TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg
end

function computeEmissionAllExplicit()

    solPlanner = optimPlannerExplicit(SOpG,quiet=true)
    GComputed = computeG(SOpG,solPlanner)
    G1Planner = computeG1(SOpG,solPlanner)
    G2Planner = computeG2(SOpG,solPlanner)
    Dp1Planner,Dp2Planner = splitInitialData(DpG,G1Planner,G2Planner)
    solODEPlanner1 = solveODE(SOpG,Dp1Planner,G1Planner)
    solODEPlanner2 = solveODE(SOpG,Dp2Planner,G2Planner)
    pPlanner1, tPlanner1 = solODEPlanner1[1,:], solODEPlanner1[2,:]
    pPlanner2, tPlanner2 = solODEPlanner2[1,:], solODEPlanner2[2,:]
    pPlanner = pPlanner1 + pPlanner2
    tPlanner = tPlanner1 + tPlanner2

    solPlannerNoResources = optimPlannerNoResourcesExplicit(SOpG,quiet=true)
    GComputed = computeG(SOpG,solPlannerNoResources)
    G1PlannerNoResources = computeG1(SOpG,solPlannerNoResources)
    G2PlannerNoResources = computeG2(SOpG,solPlannerNoResources)
    Dp1PlannerNoResources,Dp2PlannerNoResources = splitInitialData(DpG,G1PlannerNoResources,G2PlannerNoResources)
    solODEPlannerNoResources1 = solveODE(SOpG,Dp1PlannerNoResources,G1PlannerNoResources)
    solODEPlannerNoResources2 = solveODE(SOpG,Dp2PlannerNoResources,G2PlannerNoResources)
    pPlannerNoResources1, tPlannerNoResources1 = solODEPlannerNoResources1[1,:], solODEPlannerNoResources1[2,:]
    pPlannerNoResources2, tPlannerNoResources2 = solODEPlannerNoResources2[1,:], solODEPlannerNoResources2[2,:]
    pPlannerNoResources = pPlannerNoResources1 + pPlannerNoResources2
    tPlannerNoResources = tPlannerNoResources1 + tPlannerNoResources2

    solNash = optimNashExplicit(SOpG,quiet=true)
    GComputed = computeG(SOpG,solNash)
    G1Nash = computeG1(SOpG,solNash)
    G2Nash = computeG2(SOpG,solNash)
    Dp1Nash,Dp2Nash = splitInitialData(DpG,G1Nash,G2Nash)
    solODENash1 = solveODE(SOpG,Dp1Nash,G1Nash)
    solODENash2 = solveODE(SOpG,Dp2Nash,G2Nash)
    pNash1, tNash1 = solODENash1[1,:], solODENash1[2,:]
    pNash2, tNash2 = solODENash2[1,:], solODENash2[2,:]
    pNash = pNash1 + pNash2
    tNash = tNash1 + tNash2

    solStackelberg = optimStackelbergExplicit(SOpG,quiet=true)
    GComputed = computeG(SOpG,solStackelberg)
    G1Stackelberg = computeG1(SOpG,solStackelberg)
    G2Stackelberg = computeG2(SOpG,solStackelberg)
    Dp1Stackelberg,Dp2Stackelberg = splitInitialData(DpG,G1Stackelberg,G2Stackelberg)
    solODEStackelberg1 = solveODE(SOpG,Dp1Stackelberg,G1Stackelberg)
    solODEStackelberg2 = solveODE(SOpG,Dp2Stackelberg,G2Stackelberg)
    pStackelberg1, tStackelberg1 = solODEStackelberg1[1,:], solODEStackelberg1[2,:]
    pStackelberg2, tStackelberg2 = solODEStackelberg2[1,:], solODEStackelberg2[2,:]
    pStackelberg = pStackelberg1 + pStackelberg2
    tStackelberg = tStackelberg1 + tStackelberg2

    return pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2

end



function computeEmissionAll()

    solPlanner = optimPlanner(SOpG,quiet=true)
    GComputed = computeG(SOpG,solPlanner)
    G1Planner = computeG1(SOpG,solPlanner)
    G2Planner = computeG2(SOpG,solPlanner)
    Dp1Planner,Dp2Planner = splitInitialData(DpG,G1Planner,G2Planner)
    solODEPlanner1 = solveODE(SOpG,Dp1Planner,G1Planner)
    solODEPlanner2 = solveODE(SOpG,Dp2Planner,G2Planner)
    pPlanner1, tPlanner1 = solODEPlanner1[1,:], solODEPlanner1[2,:]
    pPlanner2, tPlanner2 = solODEPlanner2[1,:], solODEPlanner2[2,:]
    pPlanner = pPlanner1 + pPlanner2
    tPlanner = tPlanner1 + tPlanner2

    solPlannerNoResources = optimPlannerNoResources(SOpG,quiet=true)
    GComputed = computeG(SOpG,solPlannerNoResources)
    G1PlannerNoResources = computeG1(SOpG,solPlannerNoResources)
    G2PlannerNoResources = computeG2(SOpG,solPlannerNoResources)
    Dp1PlannerNoResources,Dp2PlannerNoResources = splitInitialData(DpG,G1PlannerNoResources,G2PlannerNoResources)
    solODEPlannerNoResources1 = solveODE(SOpG,Dp1PlannerNoResources,G1PlannerNoResources)
    solODEPlannerNoResources2 = solveODE(SOpG,Dp2PlannerNoResources,G2PlannerNoResources)
    pPlannerNoResources1, tPlannerNoResources1 = solODEPlannerNoResources1[1,:], solODEPlannerNoResources1[2,:]
    pPlannerNoResources2, tPlannerNoResources2 = solODEPlannerNoResources2[1,:], solODEPlannerNoResources2[2,:]
    pPlannerNoResources = pPlannerNoResources1 + pPlannerNoResources2
    tPlannerNoResources = tPlannerNoResources1 + tPlannerNoResources2

    solNash = optimNash(SOpG,quiet=true)
    GComputed = computeG(SOpG,solNash)
    G1Nash = computeG1(SOpG,solNash)
    G2Nash = computeG2(SOpG,solNash)
    Dp1Nash,Dp2Nash = splitInitialData(DpG,G1Nash,G2Nash)
    solODENash1 = solveODE(SOpG,Dp1Nash,G1Nash)
    solODENash2 = solveODE(SOpG,Dp2Nash,G2Nash)
    pNash1, tNash1 = solODENash1[1,:], solODENash1[2,:]
    pNash2, tNash2 = solODENash2[1,:], solODENash2[2,:]
    pNash = pNash1 + pNash2
    tNash = tNash1 + tNash2

    solStackelberg = optimStackelberg(SOpG,quiet=true)
    GComputed = computeG(SOpG,solStackelberg)
    G1Stackelberg = computeG1(SOpG,solStackelberg)
    G2Stackelberg = computeG2(SOpG,solStackelberg)
    Dp1Stackelberg,Dp2Stackelberg = splitInitialData(DpG,G1Stackelberg,G2Stackelberg)
    solODEStackelberg1 = solveODE(SOpG,Dp1Stackelberg,G1Stackelberg)
    solODEStackelberg2 = solveODE(SOpG,Dp2Stackelberg,G2Stackelberg)
    pStackelberg1, tStackelberg1 = solODEStackelberg1[1,:], solODEStackelberg1[2,:]
    pStackelberg2, tStackelberg2 = solODEStackelberg2[1,:], solODEStackelberg2[2,:]
    pStackelberg = pStackelberg1 + pStackelberg2
    tStackelberg = tStackelberg1 + tStackelberg2

    return pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2

end



function plotTemperature(;explicit=true)
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    if explicit & (SOpG.α == 1)
        TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg = computeTemperatureAllExplicit()        
    else
        TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg = computeTemperatureAll()
    end

    plot(DpG.TSave, TempPlanner, label = "Global planner",color="blue")
    plot(DpG.TSave, TempPlannerNoResources, label = "Restricted planner",color="gold")
    plot(DpG.TSave, TempNash, label = "Nash",color="green")
    plot(DpG.TSave, TempStackelberg, label = "Stackelberg",color="red")
    xlabel("Time (years)")
    ylabel("Global temperature (°C)")
    grid()
    legend()
    tight_layout()

    savefig("Temperature.pdf")

end

function plotTotalEmission(;explicit = true)
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    if explicit & (SOpG.α == 1)
        pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2 = computeEmissionAllExplicit()
    else
        pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2 = computeEmissionAll()
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

    plot(DpG.TSave, pPlanner + tPlanner, label = "Global planner - total emissions",color="blue")
    plot(DpG.TSave, pPlanner1 + tPlanner1, label = "Global planner - country 1 emissions",color="blue", "--")
    plot(DpG.TSave, pPlanner2 + tPlanner2, label = "Global planner - country 2 emissions",color="blue", ":")
    plot(DpG.TSave, pPlannerNoResources + tPlannerNoResources, label = "Restricted planner - total emissions",color="gold")
    plot(DpG.TSave, pPlannerNoResources1 + tPlannerNoResources1, label = "Restricted planner - country 1 emissions",color="gold","--")
    plot(DpG.TSave, pPlannerNoResources2 + tPlannerNoResources2, label = "Restricted planner - country 2 emissions",color="gold",":")
    plot(DpG.TSave, pNash + tNash, label = "Nash - total emissions",color="green")
    plot(DpG.TSave, pNash1 + tNash1, label = "Nash - country 1 emissions",color="green","--")
    plot(DpG.TSave, pNash2 + tNash2, label = "Nash - country 2 emissions",color="green",":")
    plot(DpG.TSave, pStackelberg + tStackelberg, label = "Stackelberg - total emissions",color="red")
    plot(DpG.TSave, pStackelberg1 + tStackelberg1, label = "Stackelberg - country 1 emissions",color="red","--")
    plot(DpG.TSave, pStackelberg2 + tStackelberg2, label = "Stackelberg - country 2 emissions",color="red",":")

    xlabel("Time (years)")
    ylabel("CO2 concentration (ppm)")
    grid()
    legend()
    tight_layout()

    savefig("TotalEmissions.pdf")
end

function plotPTEmission(;explicit = true)
    fig = gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)

    if explicit & (SOpG.α == 1)
        pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2 = computeEmissionAllExplicit()
    else
        pPlanner, tPlanner, pPlanner1, tPlanner1, pPlanner2, tPlanner2, pPlannerNoResources, tPlannerNoResources, pPlannerNoResources1, tPlannerNoResources1, pPlannerNoResources2, tPlannerNoResources2, pNash, tNash, pNash1, tNash1, pNash2, tNash2, pStackelberg, tStackelberg, pStackelberg1, tStackelberg1, pStackelberg2, tStackelberg2 = computeEmissionAll()
    end

    plot(DpG.TSave, pPlanner, label = "Global planner - permanent emission",color="blue","-")
    plot(DpG.TSave, pPlannerNoResources, label = "Restricted planner - permanent emission",color="gold","-")
    plot(DpG.TSave, pNash, label = "Nash - permanent emission",color="green","-")
    plot(DpG.TSave, pStackelberg, label = "Stackelberg - permanent emission",color="red","-")
    legend(loc=7, title="(Left scale)")    
    grid()
    TemporaryPlot = twinx()
    TemporaryPlot.tick_params(axis="y",color="red",labelcolor="red")
    
    plot(DpG.TSave, tPlanner, label = "Global planner - temporary emission",color="blue","--")
    plot(DpG.TSave, tPlannerNoResources, label = "Restricted planner - temporary emission",color="gold","--")
    plot(DpG.TSave, tNash, label = "Nash - temporary emission",color="green","--")
    plot(DpG.TSave, tStackelberg, label = "Stackelberg - temporary emission",color="red","--")

    xlabel("Time (years)")
    ylabel("CO2 concentration (ppm)")

    legend(loc=9,title="(Right scale)") 
    tight_layout()

    savefig("PermTempEmissions.pdf")

end



function plotgh(;explicit = true)

    if explicit & (SOpG.α == 1)
        allPlanner = computeAllPlannerExplicit()
        allPlannerNoResources = computeAllPlannerNoResourcesExplicit()
        allNash = computeAllNashExplicit()
        allStackelberg = computeAllStackelbergExplicit()
    else
        allPlanner = computeAllPlanner()
        allPlannerNoResources = computeAllPlannerNoResources()
        allNash = computeAllNash()
        allStackelberg = computeAllStackelberg()
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

# G = G1 + G2

# P' = ϕ G
# P1' = ϕ G1
# P2' = ϕ G2

# T' = -ϕT + ϕL * G
# T1' = -ϕT1 + ϕL * G1
# T2' = -ϕT2 + ϕL * G2

# dividi T0 e P0 come G1 e G2

# show S1 = P1 + T1 
# show S2 = P2 + T2
# S = S1 + S2

# Temp = 3/log(2)*log(S/S2000) 
# Temp1 = S1/S * 3/log(2)*log(S/S2000) 
# Temp2 = S2/S * 3/log(2)*log(S/S2000) 

# show Temperature1 
# show Temperature2
