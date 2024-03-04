
const varNames = ["C1","C2","B1","B2","K1","K2","Ra","Rb"]
const GLOBAL_TOL = 1e-16
const GLOBAL_MAX_IT = Int(1e6)

@with_kw mutable struct climateSOParameters{T}
    # static optimization parameters

    # Production
    α::T = 1.0
    A̅::T = 6.0

    # Utility
    σ1::T = 1.0
    σ2::T = 1.0
    γ1::T = 15.0 * 1e-3 * 0.5
    γ2::T = 35.0 * 1e-3 * 0.5
    # γ1::T = 20.0 * 1e-3 * 0.5
    # γ2::T = 20.0 * 1e-3 * 0.5
    ρ::T = -log((0.96)^10)

    # robustness    
    γ̂1::T = γ1      # gamma1 nature
    γ̂2::T = γ2      # gamma2 nature
    αR::T = 1e5     # alpha nature (alpha1 = alpha2)
    # αR::T = 1e4         # alpha nature (alpha1 = alpha2)

    # Pollution
    ϕL::T = 0.2
    ϕ::T = 0.5
    ϕ0::T = 0.393

    # Depollution
    θ1::T = 0.5
    θ2::T = 0.5
    ηK::T = 1.0
    ηB::T = 1.0

    h0::T = 0.5
    h∞::T = 0.9
    g0::T = 0.2
    g∞::T = 0.5

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
    # @assert 0 < δI <= 1
    @assert gθ_isConcave(g0,g∞,θ2)

end


@with_kw mutable struct climateDynamicParameters{T}    
    TInitial::T = 2000.0        # initial year
    TSim::T = 20.0             # how many years to simulate
    TSave::LinRange{T,Int64} = LinRange(TInitial, TInitial+TSim, 1000)
    P0::T = 684.0               # permanent emission at time T0
    T0::T = 118.0               # temporary emission at time T0
    S̅::T = 581.0                # Preindustrial GHC
end

function G(K1,K2,B1,B2,Rb,p)
    return p.ηK * (K1+K2) - p.ηB * (B1^p.θ1 + g(Rb,p) * B2^p.θ2)
end

function G1(K1,B1,p)
    return p.ηK * K1 - p.ηB * ( B1^p.θ1 )
end

function G2(K2,B2,Rb,p)
    return p.ηK * K2 - p.ηB * ( g(Rb,p) * B2^p.θ2 )
end


function g(x,p)
    return (p.g∞ * x + p.g0)/(x+1)
end

function h(x,p)
    return (p.h∞ * x + p.h0)/(x+1)
end

function computeG(p,sol)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    return G(K1,K2,B1,B2,Rb,p)
end

function computeG1(p,sol)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    return G1(K1,B1,p)
end

function computeG2(p,sol)
    C1,C2,B1,B2,K1,K2,Ra,Rb = [sol[name] for name in varNames]
    return G2(K2,B2,Rb,p)
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

function print_rgb(r, g, b, t)
    print("\e[1m\e[38;2;$r;$g;$b;249m",t)
end