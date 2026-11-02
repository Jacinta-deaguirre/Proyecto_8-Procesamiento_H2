using Symbolics

function test()
    @variables T nₐ
    # Modified Dubinin-Astakov Isotherm parameters
    α = 3080.0 # Enthalpic Factor / J/mol
    β = 18.9 # Entropic Factor / J/mol K
    m = 2.0 # Exponential factor
    p₀ = 1470e6 # Saturation pressure / Pa
    n₀ = 71.6 # Limit adsoption / mol/kg
    R = 8.314 # Universal gas constant / J/(mol K)  

    # Define de symbolic_expresion
    logP(T, nₐ) = log(p₀) - (α / (R * T) + β / R) * log(n₀ / nₐ)^(1 / m)

    symbolic_expresion = logP(T, nₐ)

    dh = Symbolics.derivative(symbolic_expresion, T) * R * T^2
    println("dh: ", dh)

    ΔH = build_function(dh, [T, nₐ], expression=Val{false})
    ΔH([281, 0.009939043824520061])
end

function test_DA()
    @variables T P
    # Dubinin-Astakov parameters
    P_lim = 77.75e6     # Pa
    ψ = 7.3235          # mmol g-1
    β = -0.0088         # mol kg-1 K-1
    κ = 772.92          # J mol-1
    γ = 18.828
    m = 2.0            # Exponent in the isotherm equation
    R = 8.314          # J/(mol K)

    nₐ(T, P) = (ψ + β * T) * exp((R * T / (κ + γ * T))^m * log(P_lim / P)^m)
    symbolic_expresion = nₐ(T, P)

    dnₐ_dT = Symbolics.derivative(symbolic_expresion, T)
    dnₐ_dP = Symbolics.derivative(symbolic_expresion, P)

    dH = -R * T / P * dnₐ_dT / dnₐ_dP
    ΔH = build_function(dH, [T, P], expression=Val{false})

    ΔH([295, 0.102564103e6])
end

test_DA()