# Started on April 26th 2023
# Cristiano Ricci - cristiano.ricci6@gmail.com

include("L_LoadAll.jl")

SOpG = climateSOParameters()
DpG = climateDynamicParameters()

function testSO()

    # planner
    optimPlanner(SOpG,quiet=true)
    optimPlannerExplicit(SOpG,quiet=true)
    optimPlannerRobust(SOpG,DpG,quiet=true)[1]

    # restricted planner
    optimPlannerNoResources(SOpG,quiet=true)
    optimPlannerNoResourcesExplicit(SOpG,quiet=true)
    optimPlannerNoResourcesRobust(SOpG,DpG,quiet=true)[1]

    # nash
    optimNash(SOpG,quiet=true)
    optimNashExplicit(SOpG,quiet=true)
    optimNashRobust(SOpG,DpG,quiet=true)[1]

    # stackelberg
    optimStackelberg(SOpG,quiet=true)
    optimStackelbergExplicit(SOpG,quiet=true)
    optimStackelbergRobust(SOpG,DpG,quiet=true)[1]

    return nothing
end

function testG()

    @show computeGPlanner(SOpG)
    @show computeGPlannerNoResources(SOpG)
    @show computeGNash(SOpG)
    @show computeGStackelberg(SOpG)
    return nothing
end

