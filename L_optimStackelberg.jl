function optimStackelberg(p; quiet = false)
    # initial point by nash optimum
    solPlanner = optimNash(p; quiet = true)
    B1,I1,Ra,Rb = solPlanner["B1"],solPlanner["I1"],solPlanner["Ra"],solPlanner["Rb"]
    
    C2,B2,I2 = Br2(p,B1,I1,Ra,Rb)
    C1,B1,I1,Ra,Rb = interate1(p,C2,B2,I2)
    candidateSol = [C1,B1,I1,Ra,Rb]

    err = 1.0; Nit = 0
    while (err > TOL) & (Nit < MAX_IT)
        
        C2,B2,I2 = Br2(p,B1,I1,Ra,Rb)
        C1,B1,I1,Ra,Rb = iterate1(p,C2,B2,I2)
        
        err = norm(candidateSol - [C1,B1,I1,Ra,Rb],Inf)
        Nit = Nit + 1 

        candidateSol = [C1,B1,I1,Ra,Rb]
        
        if quiet == false @show err, Nit end
    end

    C2,B2,I2 = Br2(p,B1,I1,Ra,Rb)
    solValues = [C1,C2,B1,B2,I1,I2,Ra,Rb]    

    varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]
    solDict = Dict(zip(varNames,solValues))
    return solDict
end
    
function Br2(p,B1,I1,Ra,Rb)

    return C2,B2,I2
end

function iterate1(p,C2,B2,I2)

    return C1,B1,I1,Ra,Rb    
end