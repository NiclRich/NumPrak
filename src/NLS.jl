using SparseArrays
using LinearAlgebra

"""
    dfp(initial_values::Function, T::Float64, dt::Float64, dx::Float64; lambda::Float64=2, a::Float64=-30, b::Float64=30)

Solves the 1D Nonlinear Schrödinger equation using a finite difference scheme.

# Arguments
- `initial_values`: A function that generates the initial values based on spatial points.
- `T`: The total simulation time.
- `dt`: The time step size.
- `dx`: The spatial step size (discretization in space).
- `lambda`: The nonlinearity parameter (default: `2`).
- `a`: The left boundary of the spatial domain (default: `-30`).
- `b`: The right boundary of the spatial domain (default: `30`).

# Returns
- `U`: A matrix where each column represents the solution at a time step.
	To be more precise U[:, k] contains the values at time (k-1) * dt

# Notes
- The time interval is [0, T].
- The spatial domain is [a, b].
- This scheme solves the Equation 

\$\$
i u_t - u'' +  |u|^2u = 0
u(x,0) = \\phi(x)
\$\$

with Dirichlet boundary conditions `` u(a,t) = 0 = u(b,t)`` for all \$ t\$  .

- The function implements a fixed-point iteration to solve for the solution at each time step.	The paper [1] does not prove the convergence of the scheme nor of the fixed-point iteration.

# References
[1] M. Delfour, M. Fortin, G. Payre: Finite Difference Solutions of a Non-linear Schrödinger
equation. Journal of Computational Physics 44, 277-288 (1981)

"""
function dfp(initial_values::Function, T, dt, dx; lambda=-2, a=-30, b=30)
	max_iteration = 1000 	# maximum number of iterations in fix point loop
	eps = 1e-14		# break criteria for fix point loop
	
	# compute the number points  in space such that (b-a) / dx is an 
	# integer. The stepsize dx will be corrected and sightly decreased.
	J = ceil(Int64, 0.5 * (b-a) / dx) # compute the number of cells
	dx = (b-a) / J 

	# Compute the number of points in the time dimension. The implicit
	# assumption is that the time starts at 0
	K = ceil(Int64, T / dt)
	dt = T / K 

	# DEBUGGING:
	# J = Int(30 / dx)
	# K = Int(6 / dt)

	# size of the quadratic sparse matrix 
	matrix_size = 2 * J -1

	# create the sparse matrix A
	A = spdiagm(0  => 2 * ones(matrix_size),  # main diagonal
		    -1 => -1 * ones(matrix_size - 1), # lower minor diagonal
		    1  => -1 * ones(matrix_size -1 ))  # upper minor diagonal
	# create sparse identity matrix
	I = sparse(1:matrix_size,1:matrix_size,ones(matrix_size))
	
	# create non linearity approximation
	F(U, Uk, lambda) = (lambda / 4) * (abs2.(U) + abs2.(Uk)) .* (U + Uk)

	# Define iteration matrices.
	# I is the identity matrix with approbiate size.
	# Mp1 for the matrix in front of U^{k+1}
	# Mk for the matrix in fron of U^k
	# Mp1 = spdiagm(0 => im / dt * ones(matrix_size)) .+ (1 / (2 * dx^2)) * A
	# Mk = spdiagm(0 => im / dt * ones(matrix_size)) .- (1 / (2 * dx^2)) * A
	Mplus = im / dt * I .+ (1 / (2 * dx^2)) * A
	Mminus = im / dt * I .- (1 / (2 * dx^2)) * A

	# Factorization for efficient solving
	FMplus = factorize(Mplus)

	# create initial values
	# x = collect(range(-30, stop=30, length=matrix_size))
	# U0 = initial_values.(x)
	U = Array{ComplexF64}(undef, 2*J-1, K+1)

    #initialise the initial condition
    for j in eachindex(U[:,1])
        U[j, 1] = initial_values((j-J)*dx)
    end

	# Preallocate the solution matrix for performance
	# U = zeros(ComplexF64, matrix_size, K+1)
	# U[:, 1] = U0

	# main loop for computing the values
	# outer loop is the time advancement step
	for k in 1:K
		error = 1 + eps
		Up = U[:, k]
		iterations = 0
		while error > eps
			Up1 = FMplus \ (Mminus * U[:, k] .- F(Up, U[:, k], lambda))
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

# exact_sol(x,t) = 1.5 * exp(im * (2 * x - 1.75 * t)) * sech(1.5 * (x+5) - 6*t)
# init_val(x) = exact_sol(x, 0)
# T = 6
# dt = 0.02
# dx = 0.1
# dxs = [10.0^(-k) for k in 0:4]
# dts = [10.0^(-k) for k in 0:4]
#  U = dfp(init_val, T,dt, dx)
#  xs = collect(range(start=-30, stop=30, length=size(U,1)))
#  ts = collect(range(0, T, length=size(U,2)))
#  xs = reshape(xs, length(xs), 1)
#  ts = reshape(ts, 1, length(ts))

#  U_true = exact_sol.(xs, ts)
#  error = maximum(abs.(U - U_true))
#  print(U)
#  @info error

#  @info "done"

