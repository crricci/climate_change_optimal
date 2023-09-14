
const varNames = ["C1","C2","B1","B2","I1","I2","Ra","Rb"]

@with_kw mutable struct climateSOParameters{T}
    # static optimization parameters

    # Production
    α::T = 2/3
    A̅::T = 1000.0

    # Utility
    σ1::T = 2.0
    σ2::T = 2.0
    γ1::T = 25.0 
    γ2::T = 25.0
    # γ1::T = 1e-2
    # γ2::T = 1e-2
    ρ::T = 0.04

    # Pollution
    ϕL::T = 0.2
    ϕ::T = 0.5
    ϕ0::T = 0.393

    # Depollution
    θ1::T = 0.5
    θ2::T = 0.5
    cI::T = 10.0
    cD::T = 0.001

    h0::T = 1/3
    h∞::T = 0.9
    g0::T = 0.5
    g∞::T = 0.9

    b̅::T = 1.0      # fixed t 1
    δI::T = 1.0     # fixed to 1
    δD::T = 1.0     # fixed to 1

    # Derived
    Φ::T = (ϕL / ρ + (1-ϕL) * ϕ0 / (ρ + ϕ))
    gPrime0::T = g∞ - g0

    # CHECKS
    @assert 0 < α <= 1
    @assert A̅*h0 > 1
    @assert h0 < 1
    @assert h∞ < 1
    @assert h0 < h∞
    @assert g0 < 1
    @assert g∞ < 1
    @assert g0 < g∞
    @assert θ1 < 1
    @assert θ2 < 1
    @assert 0 < δI <= 1
    @assert δD >= 1
    @assert gθ_isConcave(g0,g∞,θ2)

end


@with_kw mutable struct climateDynamicParameters{T}    
    TInitial::T = 2000.0        # initial year
    TSim::T = 100.0             # how mazny years to simulate
    TSave::LinRange{T,Int64} = LinRange(TInitial, TInitial+TSim, 1000)
    P0::T = 684.0               # permanent emission at time T0
    T0::T = 118.0               # temporary emission at time T0
    S̅::T = 581.0                # Preindustrial GHC
end

function G(I1,I2,B1,B2,Rb,p)
    return p.cI * (I1+I2)^p.δI - p.cD * (p.b̅ * B1^p.θ1 + p.b̅ * g(Rb,p) * B2^p.θ2)^p.δD
end

function g(x,p)
    return (p.g∞ * x + p.g0)/(x+1)
end

function h(x,p)
    return (p.h∞ * x + p.h0)/(x+1)
end

function computeG(p,sol)

    C1,C2,B1,B2,I1,I2,Ra,Rb = [sol[name] for name in varNames]
    return G(I1,I2,B1,B2,Rb,p)
end

function gθ_isConcave(g0,g∞,θ)

    x = LinRange(0,100,1000);

    Symbolics.@variables xS g0S g∞S θS;
    fExpr = ((g∞S * xS + g0S)/(xS+1))^(1/(1-θS));
    D = Differential(xS);
    fExpr2nd = expand_derivatives(D(D(fExpr)))

    if any(p -> p >= 0, [substitute(fExpr2nd,(Dict(xS => xi, g0S => g0, g∞S => g∞, θS => θ))) for xi in x])
        return false
    else
        return true
    end
end

