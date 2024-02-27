
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
    return 3/log(2) * log((P+T)/Dp.S̅)
end

function splitInitialData(Dp,G1,G2)

    Dp1 = deepcopy(DpG)
    Dp2 = deepcopy(DpG)

    Dp1.P0 = Dp.P0 / 2
    Dp2.P0 = Dp.P0 / 2

    Dp1.T0 = Dp.T0 / 2
    Dp2.T0 = Dp.T0 / 2
    

    # if (G1 > 0) & (G2 > 0 )
    #     Dp1.P0 = Dp.P0 * G1 / (G1+G2) 
    #     Dp1.T0 = Dp.T0 * G1 / (G1+G2) 
    #     Dp2.P0 = Dp.P0 * G2 / (G1+G2) 
    #     Dp2.T0 = Dp.T0 * G2 / (G1+G2) 
    # elseif (G1 > 0) & (G2 < 0)
    #     Dp2.P0 = 0.0
    #     Dp2.T0 = 0.0
    # elseif (G1 < 0) & (G2 > 0)
    #     Dp1.P0 = 0.0
    #     Dp1.T0 = 0.0
    # else
    #     error("G total cannot be negative")
    # end

    return Dp1,Dp2 
    
end