# Started on April 26th 2023
# Cristiano Ricci - cristiano.ricci6@gmail.com

include("L_LoadAll.jl")

function main(p)

    @show optimPlanner(p,quiet=true)
    @show optimPlannerNoResources(p,quiet=true)
    @show optimNash(p,quiet=true)

end

