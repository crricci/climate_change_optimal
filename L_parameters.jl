
@with_kw mutable struct climateFramework{T}

    # Production
    α::T = 1.0
    A̅::T = 3.0
    h0::T = 0.5
    h∞::T = 0.9

    # Utility
    σ1::T = 1.0
    σ2::T = 1.0
    k::T = 1.0
    γ1::T = 1.0
    γ2::T = 1.0
    ρ::T = 0.5

    # Pollution
    ϕL::T = 1.0
    ϕ::T = 1.0
    ϕ0::T = 1.0
    P0::T = 1.0
    T0::T = 1.0

    # Depollution
    b̅::T = 1.0
    g0::T = 0.2
    g∞::T = 0.9
    θ1::T = 0.5
    θ2::T = 0.5
    cI::T = 1.0
    cD::T = 1.0
    δI::T = 1.0
    δD::T = 1.0

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

