using LinearAlgebra
using SparseArrays
using Plots

function schroedingereq(init_val::Function, dx, dt; a=-30, b=30)
    epsilon = 1e-7
    lambda = -2
    T = 6.0

    # # compute the number points  in space such that (b-a) / dx is an 
	# # integer. The stepsize dx will be corrected and sightly decreased.
	# J = ceil(Int64, 0.5 * (b-a) / dx) # compute the number of cells
	# dx = (b-a) / J 

	# # Compute the number of points in the time dimension. The implicit
	# # assumption is that the time starts at 0
	# K = ceil(Int64, T / dt)
	# dt = T / K 

    J = Int(30.0 / dx)
    T  = 6.0
    K = Int(T / dt)

    # define matrices

    m_size = 2 * J - 1
    A = spdiagm(0 => 2 * ones(m_size),
                -1 => -1 * ones(m_size - 1),
                1 => -1 * ones(m_size - 1))

    I = sparse(1:m_size, 1:m_size, ones(m_size))

    Mplus = im / dt * I + 1 / (2 * dx * dx) * A
    Mminus = im / dt * I - 1 / (2 * dx * dx) * A

    # create non linearity approximation
	F(U, Uk, lambda) = (lambda / 4) * (abs2.(U) + abs2.(Uk)) .* (U + Uk)

    # define solution matrix
    U = zeros(ComplexF64, m_size, K+1)

    xs = collect(range(start=-30, stop=30, length=m_size))
    U[:, 1] = init_val.(xs)

    # main loop for computing the values
	# outer loop is the time advancement step
    max_iteration = 1000
	for k in 1:K
		error = 1 + epsilon
		Up = deepcopy(U[:, k])
        Uk = deepcopy(U[:, k])
		iterations = 0
		while error > epsilon
            Fk = 0.25 * lambda * (abs2.(Uk) + abs2.(Up)) .*(Uk + Up)
			Up1 = Mplus \ (Mminus * Uk - Fk)
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
    @info "Computation done"
    return U
end
# exact_sol(x,t) = 1.5 .* exp.(im .* (2.0 * x .- 1.75 .* t)) .* sech.(1.5 .* (x .+ 5) .+ 6 .* t)

exact_sol(x,t) = 1.5 .* exp.(im .* (1.25 .* t - x )) .* sech.(1.5 .* (x .+ 5) .- 3.0 * t)
init(x) = exact_sol(x, 0)
# dxs = [1, 0.5, 0.1, 0.05, 0.01]
# dts = [1, 0.5, 0.1, 0.05, 0.01]
# for dt in dts

dt = 0.01
dx = 0.1
U_num = schroedingereq(init, dx, dt)

ts = collect(range(start=0, stop=6.0, length=Int(6.0 / dt)+1))
xs = collect(range(start=-30, stop=30.0, length=Int(30 / dx) * 2 - 1))
xs = reshape(xs, length(xs), 1)
ts = reshape(ts, 1, length(ts))
U_true = exact_sol.(xs, ts)
@info size(exact_sol.(xs, ts))
@info size(U_num)

@info "Fehler: " norm(U_true - U_num, Inf)
@info "Fehler abs" abs(norm(U_true, Inf) - norm(U_num, Inf))
# end

# plot([norm(abs.(U_num[:, k] - U_true[:, k]), Inf) for k in 1:300])
# plot(xs, [imag.(U_num[:, 50]), abs.(U_true[:, 50])])
# @gif for k in 1:300
#     plot(xs, [abs.(U_num[:, k]), abs.(U_true[:, k])])
# end

# Check the conservation laws
first_conservation = sum(abs2, U_num, dims=1)
out = sum([1 2 3])
plot(ts, first_conservation)
