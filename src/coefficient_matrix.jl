"""
Function that computes the coefficient matrix for a radial heat conduction component.
It starts from the component of heat conduction in cylindrical coordinates:

    k/r ∂/∂r(r ∂T/∂r)

By using the finite difference method, we can discretize this equation and by rearranging
we can obtain the coefficient matrix for the system of equations that describes the heat conduction in a radial direction:

    AT + b

Inputs:
- `R`: Radius of the cylinder (m)
- `dr`: Radial step size (m)
- `k_eff`: Effective thermal conductivity (W/(m·K))
- `U`: Heat transfer coefficient (W/(m²·K))
- `T_infty`: Ambient temperature (K)

Outputs:
- `n_r`: Number of radial nodes
- `r`: Radial positions vector of dimension n_r
- `A`: Coefficient matrix of dimension n_r x n_r 
- `b`: Right-hand side vector of dimension n_r
"""
function coefficient_matrix(R::Real, dr::Real, k_eff::Real, U::Real, T_infty::Real)
    n_r = floor(Int32, R / dr) + 1 # Number of radial nodes
    r = range(0, stop=R, length=n_r) # Radial positions

    # Initialize the coefficient matrix and the right-hand side vector
    A = zeros(n_r, n_r)
    b = zeros(n_r)

    # Constants for the finite difference method
    k1 = k_eff / (2 * dr)
    k2 = k_eff / (dr^2)

    # fill the matrix internal nodes  
    for i in 2:n_r
        #println("i: $i, r: $((i - 1) * dr)")
        if i < n_r
            A[i, i-1] = -k1 / r[i] + k2
            A[i, i] = -2 * k2
            A[i, i+1] = k1 / r[i] + k2
        end
    end

    # Boundary conditions
    # Axial symmetry
    #A[1, 1] = -4 * k2
    #A[1, 2] = 4 * k2
    A[1, 1] = -3
    A[1, 2] = 4
    A[1, 3] = -1

    # Heat transfer to the exterior
    #A[n_r, n_r] = 3 + 2 * U * dr / k_eff
    #A[n_r, n_r-1] = -4
    #A[n_r, n_r-2] = 1

    # b[n_r] = 2 * dr * U / k_eff * T_infty # Should be 0 because it's been moved to T_wall
    return n_r, r, A, b
end