"""
Function that computes the coefficient matrix for a radial heat conduction component.
It starts from the component of heat conduction in cylindrical coordinates:

    k/r ∂/∂r(r ∂T/∂r)

By using the finite difference method, we can discretize this equation and by rearranging
we can obtain the coefficient matrix for the system of equations that describes the heat conduction in a radial direction:

    AT + b

Inputs:
- `R`: Radius of the cylinder (m)
- `dr`: Nominal radial step size (m). The number of nodes is `n_r = round(R/dr) + 1` and the
        grid `r = range(0, R, length = n_r)` covers the radius exactly; the actual spacing
        `step(r)` (which differs from `dr` by less than 0.1 %) is the one used in the
        finite-difference coefficients. Scripts must use `r[2] - r[1]` as the spacing.
- `k_eff`: Effective thermal conductivity (W/(m·K))
- `U`: Heat transfer coefficient (W/(m²·K)) (kept for backward compatibility, unused)
- `T_infty`: Ambient temperature (K) (kept for backward compatibility, unused)

Outputs:
- `n_r`: Number of radial nodes
- `r`: Radial positions vector of dimension n_r
- `A`: Coefficient matrix of dimension n_r x n_r
- `b`: Right-hand side vector of dimension n_r
"""
function coefficient_matrix(R::Real, dr::Real, k_eff::Real, U::Real, T_infty::Real)
    n_r = round(Int, R / dr) + 1 # Number of radial nodes
    r = range(0, stop=R, length=n_r) # Radial positions (exact radius)
    dr_a = step(r) # Actual radial spacing used in the finite differences

    # Initialize the coefficient matrix and the right-hand side vector
    A = zeros(n_r, n_r)
    b = zeros(n_r)

    # Constants for the finite difference method (second-order central differences)
    k1 = k_eff / (2 * dr_a)
    k2 = k_eff / (dr_a^2)

    # fill the matrix internal nodes
    for i in 2:n_r-1
        A[i, i-1] = -k1 / r[i] + k2
        A[i, i] = -2 * k2
        A[i, i+1] = k1 / r[i] + k2
    end

    # Boundary conditions
    # Axial symmetry at r = 0: second-order one-sided difference of ∂T/∂r = 0,
    # -3 T_1 + 4 T_2 - T_3 = 0. The diffusion operator is not evaluated at the axis;
    # the row is used through its time derivative in the DAE residual.
    A[1, 1] = -3
    A[1, 2] = 4
    A[1, 3] = -1

    # Bed-wall interface (r = R): handled as an algebraic equation in the DAE residual
    # (flux continuity with the lumped wall), so the last row is left empty here.
    return n_r, r, A, b
end
