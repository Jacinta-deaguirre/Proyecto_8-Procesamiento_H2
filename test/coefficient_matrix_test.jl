# Test parameters
R = 0.5  # Radius of the cylinder (m)
dr = 0.1  # Radial step size (m)
k_eff = 1.0  # Effective thermal conductivity (W/(m·K))
U = 10.0  # Heat transfer coefficient (W/(m²·K))
T_0 = 300.0  # Initial temperature of the tank (K)
T_infty = 350.0  # Ambient temperature (K)

# Call the function to compute the coefficient matrix
A, T₀, n_r, b, r = coefficient_matrix(R, dr, k_eff, U, T_0, T_infty)

println("n_r: ", n_r)
display(A)
display(Tᵢ)
println("b: ", b)
println("r_span: ", r_span)