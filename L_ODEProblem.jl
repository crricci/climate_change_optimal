
function solveODE(SOp,Dp,GComputed)
    tspan = (Dp.TInitial, Dp.TInitial+Dp.TSim)
    u0 = [Dp.P0, Dp.T0]

    prob = ODEProblem(df!,u0,tspan,(SOp,Dp,GComputed))
    sol = DifferentialEquations.solve(prob,Tsit5(),saveat=Dp.TSave)

    return sol
end

function df!(du,u,p,t)

    SOp,Dp,GComputed = p
    P,T = u
    dP,dT = du

    dP = SOp.ϕL * GComputed
    dT = - SOp.ϕ * T + (1- SOp.ϕL) * SOp.ϕ0 * GComputed

    du .= [dP,dT]
end

function ComputeTemperature(Dp, P,T)
    return 3*log(2) * log((P+T)/Dp.S̅)
end


function solveModelTemp()

    computeGPlanner(SOpG)
    # computeGPlannerNoResources(SOpG)
    # computeGNash(SOpG)
    # computeGStackelberg(SOpG)

    P,T = solveODE(SOpG,DpG,GComputed)
    TempFinal = ComputeTemperature(DpG,P,T)
    TempInitial = ComputeTemperature(DpG,DpG.P0,DpG.T0)
    ΔTemp = TempFinal - TempInitial
end

