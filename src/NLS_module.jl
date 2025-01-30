module NLS

export schroedingereq, explicit_NLS, first_conservation, second_conservation

using LinearAlgebra
using SparseArrays
using Plots

"""
    schroedingereq(init_val::Function, dx::Real, dt::Real; lambda::Real = -2, epsilon::Real = 1e-7, T::Real = 6.0) -> Matrix{ComplexF64}

Solve the **non-linear Schrödinger equation** using the finite-difference method described in 
*"Finite-Difference Solutions of a Non-linear Schrödinger Equation"* by Delfour, Fortin, and Payre (1981) [1].

This function implements a **Crank-Nicolson time-stepping scheme** combined with a fixed-point iteration 
to handle the non-linearity in the Schrödinger equation.

# Mathematical Context

The Schrödinger equation with non-linearity can be written as:

```math
i u_t + u_{xx} + \\lambda |u|^2 u = 0
```
with 
- `u` as the function to solve
- `i` the imaginary unit with `i^2 = -1`
- `x` as the spatial variable
- `t` as the time variable
- `\\lambda` as an equation parameter

The function discretizes the spatial domain with a step size dx and the time 
    domain with a step size dt. It uses a fixed-point iteration to approximate
    the non-linear term at each time step.

The spatial domain is fixed with the interval (-30, 30).

# Arguments
- `init_val : Function` : Function which computes the initial values
- `dx : Real` : mesh size in the spatial variable
- `dt : Real` : time step size

# Keyword Arguments
- `lambda : Real` Parameter of the NLS-Equation, default = -2.0
- `epsilon : Real` Error threshold for the fix point iteration, default = 1e-7
- `T : Real` Total simulation time, default = 6.0

# Returns
`U : Matrix{ComplexF64}` : The rows correspond to the spatial value and the columns correspond to the time values

# References# References
[1] M. Delfour, M. Fortin, G. Payre: Finite Difference Solutions of a Non-linear Schrödinger
equation. Journal of Computational Physics 44, 277-288 (1981)
"""
function schroedingereq(init_val::Function, dx, dt; lambda = -2, epsilon=1e-7, T=6.0)
    # compute the number of dimensions for the matrices
    J = Int(30.0 / dx)
    K = Int(T / dt)

    # define the size of the matrix
    m_size = 2 * J - 1

    # define the matrices for the iteration
    A = spdiagm(0 => 2 * ones(m_size),
                -1 => -1 * ones(m_size - 1),
                1 => -1 * ones(m_size - 1))

    I = sparse(1:m_size, 1:m_size, ones(m_size))
    Mplus = im / dt * I + 1 / (2 * dx * dx) * A
    Mminus = im / dt * I - 1 / (2 * dx * dx) * A

    # define solution matrix
    U = zeros(ComplexF64, m_size, K+1)

    # discretisation of the spatial variable
    xs = collect(range(start=-30, stop=30, length=m_size))

    # compute the initial values along the spatial domain
    U[:, 1] = init_val.(xs)

    # main loop for computing the values
	# outer loop is the time advancement step
    # define the maximum number of iterations
    max_iteration = 1000
	for k in 1:K
		error = Inf # to make sure, that the while loop makes atleast one iteration
		Up = deepcopy(U[:, k])
        Uk = deepcopy(U[:, k])
		iterations = 0

        # Fix Point iteration for solving the nonlinear equation
		while error > epsilon
            # compute the nonlinear term
            Fk = 0.25 * lambda * (abs2.(Uk) + abs2.(Up)) .*(Uk + Up)

            # compute the next iteration step
			Up1 = Mplus \ (Mminus * Uk - Fk)

            # compute the norm of the difference between two fix point iterations
			error = norm(Up1 .- Up)
			Up = Up1

			iterations += 1
			if iterations >= max_iteration
				@warn "Fix point iteration did not converge"
				break
			end
		end
		# Update the matrix with the solution
		U[:, k+1] = Up
    end
    return U
end


"""
    first_conservation(U::AbstractMatrix{<:Complex})

Compute the first conserved quantity of the NLS-equation.

# Description
This function calculates the **first integral**, the first integral is given by
```math
\\int |u|^2 \\; dx
```
    
It assumes:
- Rows of `U` represent the spatial variable.
- Columns of `U` represent the temporal variable.

The function computes:
1. The **first integral** as the sum of squared magnitudes (norm) of the solution matrix `U` for each time step.
2. The **maximum** and **minimum** values of the integral across the spatial domain.
3. The **dissipation**, defined as the difference between the maximum and minimum values.

# Arguments
- `U::AbstractMatrix{<:Complex}`: A complex-valued matrix where rows represent spatial variables and columns represent temporal variables.

# Returns
A tuple with the following values:
- `first_integral::Vector{Float64}`: The computed first integral for each time step.
- `maximal_value::Float64`: The maximum value of the first integral.
- `minimal_value::Float64`: The minimum value of the first integral.
- `dissipation::Float64`: The difference between the maximum and minimum values (indicating dissipation).

# Example
```julia
# Example usage
U = rand(ComplexF64, 100, 10)  # Random complex-valued matrix with 100 spatial points and 10 time steps

first_integral, max_val, min_val, dissipation = first_conservation(U)

println("First integral: ", first_integral)
println("Maximum value: ", max_val)
println("Minimum value: ", min_val)
println("Dissipation: ", dissipation)
```
"""
function first_conservation(U::Matrix{ComplexF64})
    # assumes that the rows contain the spatial variable and
    # the columns the temporal variable
    first_integral = first_conservation = sum(abs2, U, dims=1)
    maximal_value = maximum(first_integral)
    minimal_value = minimum(first_integral)
    dissipation = maximal_value - minimal_value

    return first_integral, maximal_value, minimal_value, dissipation
end

"""
    second_conservation(U::AbstractMatrix{<:Complex}, dx::Real; lambda::Real = -2)

Compute the second conserved quantity for a numerical solution of the NLS-equation.
The second conserved quantity is given by:
```math
\\int \\frac{1}{2} |u'|^2 + \\frac{\\lambda}{4}|u|^4 \\; dx
```
# Description
This function calculates the **second integral**, which typically corresponds to the energy-like conservation quantity in numerical schemes, assuming that:
- Rows of `U` represent the spatial variable.
- Columns of `U` represent the temporal variable.

The function ensures the spatial step size `dx` is positive. It computes:
1. The **second integral** as a discrete analogue of the integral.
2. The **maximum** and **minimum** values of the integral across the temporal domain.
3. The **dissipation**, defined as the difference between the maximum and minimum values.

# Arguments
- `U::AbstractMatrix{<:Complex}`: A complex-valued matrix where rows represent spatial variables and columns represent temporal variables.
- `dx::Real`: The spatial step size. Must be a positive real number.
- `lambda::Real` : Parameter of the NLS


# Returns
A tuple with the following values:
- `second_integral::Vector{Float64}`: The computed second integral for each time step.
- `maximal_value::Float64`: The maximum value of the second integral.
- `minimal_value::Float64`: The minimum value of the second integral.
- `dissipation::Float64`: The difference between the maximum and minimum values (indicating dissipation).

# Errors
- Throws an `ErrorException` if `dx` is not a positive real number.

# Example
```julia
# Example usage
U = rand(ComplexF64, 100, 10)  # Random complex-valued matrix with 100 spatial points and 10 time steps
dx = 0.1                      # Spatial step size

second_integral, max_val, min_val, dissipation = second_conservation(U, dx)

println("Second integral: ", second_integral)
println("Maximum value: ", max_val)
println("Minimum value: ", min_val)
println("Dissipation: ", dissipation)
```
"""
function second_conservation(U::Matrix{ComplexF64}, dx::Real, lambda::Real)
    # assumes that the rows contain the spatial variable and
    # the columns the temporal variable
    if dx <= 0
        error("dx must be a positive real number")
    end
    abs4(x) = abs2(x)^2


    second_integral =  1 / (2 * dx) * sum(abs2, diff(U, dims=1), dims=1) + 0.25 * lambda * dx * sum(abs4, U, dims=1)
    maximal_value = maximum(second_integral)
    minimal_value = minimum(second_integral)
    dissipation = maximal_value - minimal_value

    return second_integral, maximal_value, minimal_value, dissipation
end

"""
    explicit_NLS(init_val::Function, dx::Real, dt::Real; lambda::Real=-2, T::Real=6.0)

Solve the Nonlinear Schrödinger (NLS) equation using an explicit finite difference scheme.

# Arguments
- `init_val::Function`: A function specifying the initial condition as a function of the spatial variable.
- `dx::Real`: The spatial step size for discretization.
- `dt::Real`: The temporal step size for discretization.
- `lambda::Real` (optional): Nonlinearity parameter (default: `-2`).
- `T::Real` (optional): Total simulation time (default: `6.0`).

# Keyword Arguments
- `lambda::Float64`: Controls the strength and type of the nonlinearity in the NLS equation.
- `T`: Specifies the time duration over which the simulation runs.

# Returns
`Matrix{ComplexF64}`: A matrix where each column represents the solution vector at a specific saved time step.

# Details
The NLS equation is solved on a spatial domain `[-30, 30]` discretized with a spatial step `dx`. 
The temporal evolution is calculated with time step `dt` using an explicit finite difference scheme.

Key features include:
- Automatic adjustment of the number of saved time steps via a `proportion` parameter, enabling memory-efficient storage.
- Discretization of the spatial domain with a finite difference matrix for the second derivative.
- Handling of the nonlinear term in the NLS equation.

# Example
```julia
using LinearAlgebra, SparseArrays

# Define an initial condition
init_val(x) = exp(-x^2)

# Set parameters
dx = 0.1
dt = 0.01

# Solve the NLS equation
solution = explicit_NLS(init_val, dx, dt; lambda=-2.0, T=5.0)
```

# Notes
 - The solution matrix `U` saves the state at regular intervals determined by the proportion parameter (default: every 100th time step).
 - If a NaN is detected during the computation, a warning is issued.
"""

function explicit_NLS(init_val::Function, dx::Real, dt::Real; lambda::Real = -2.0,T::Real=6.0)

    # Only save the proportion-th number of iterations, 
    # this allows smaller time steps without running out of memory,
    # so only 1/proportion is of the time steps are saved
    proportion = 100

    # compute the number of dimensions for the matrices
    J = Int(30.0 / dx)
    K = Int(T / dt) * proportion

    # define the size of the matrix
    m_size = 2 * J - 1

    # define matrix for the 2nd derivative
    A = spdiagm(0 => 2 * ones(m_size),
                -1 => -1 * ones(m_size - 1),
                1 => -1 * ones(m_size - 1))


    # define solution matrix
    U = zeros(ComplexF64, m_size,Int(K / proportion) +1)

    # discretisation of the spatial variable
    xs = collect(range(start=-30, stop=30, length=m_size))

    # compute the initial values along the spatial domain
    U[:, 1] = init_val.(xs)
    Uk = deepcopy(U[:, 1])

    # Iteration loop for the time stepping
    for k in 1:K-1
        if k >= size(U, 2)
            break
        end
        # compute the nonlinear term
        Fk = 0.25 * lambda * (2 * abs2.(Uk)) .*(2 * Uk)
        # compute the next vector using the FD-scheme
        Uk1 = Uk + (dt / (abs2(dx) * im)) * A * Uk - (im / dt) * Fk
        if k % proportion == 0
            U[:, k+1] = Uk1
        end
        if k % proportion == 0 & k != 1
            Uk = deepcopy(U[:, k+1])
        else
            Uk = Uk1
        end
        # uncomment to detect NaN values for the unstable scheme
        # if any(isnan.(Uk1))
        #     @warn "Nan detected in iteration $k at t = $(dt * k)"
        # end
    end

    return U
end

end # end of the module
