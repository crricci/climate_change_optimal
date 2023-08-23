# Started on April 26th 2023
# Cristiano Ricci - cristiano.ricci6@gmail.com

include("L_LoadAll.jl")

function main(p)

    # planner
    @show optimPlanner(p,quiet=true)
    @show optimPlannerExplicit(p,quiet=true)

    # restricted planner
    @show optimPlannerNoResources(p,quiet=true)
    @show optimPlannerNoResourcesExplicit(p,quiet=true)

    # nash
    @show optimNash(p,quiet=true)
    @show optimNashExplicit(p,quiet=true)

    # stackelberg
    @show optimStackelberg(p,quiet=false)
    @show optimStackelbergExplicit(p,quiet=true)

end

