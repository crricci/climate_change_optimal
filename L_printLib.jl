

function makeTable()

    allStackelberg = computeAllPlanner()
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

function plotTemperature(TempPlanner, TempPlannerNoResources, TempNash, TempStackelberg)
    
    plot(DpG.TSave, TempPlanner, label = "Global planner")
    plot(DpG.TSave, TempPlannerNoResources, label = "Restricted planner")
    plot(DpG.TSave, TempNash, label = "Nash")
    plot(DpG.TSave, TempStackelberg, label = "Stackelberg")
    xlabel("Time (years)")
    ylabel("Global temperature (°C)")
    grid()
    legend()
    tight_layout()

end