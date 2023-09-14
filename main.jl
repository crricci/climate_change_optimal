# Started on April 26th 2023
# Cristiano Ricci - cristiano.ricci6@gmail.com

include("L_LoadAll.jl")

SOpG = climateSOParameters()
DpG = climateDynamicParameters()

function testSO()

    # planner
    @show optimPlanner(SOpG,quiet=true)
    @show optimPlannerExplicit(SOpG,quiet=true)

    # restricted planner
    @show optimPlannerNoResources(SOpG,quiet=true)
    @show optimPlannerNoResourcesExplicit(SOpG,quiet=true)

    # nash
    @show optimNash(SOpG,quiet=true)
    @show optimNashExplicit(SOpG,quiet=true)

    # stackelberg
    @show optimStackelberg(SOpG,quiet=true)
    @show optimStackelbergExplicit(SOpG,quiet=true)

    return nothing
end

function testG()

    @show computeGPlanner(SOpG)
    @show computeGPlannerNoResources(SOpG)
    @show computeGNash(SOpG)
    @show computeGStackelberg(SOpG)
    return nothing
end

