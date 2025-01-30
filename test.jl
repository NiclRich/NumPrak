using LinearAlgebra
using Printf
using Plots

using SparseArrays


exact_sul(x,t) = 1.5 .* exp.(im .* (1.25 .* t - x )) .* sech.(1.5 .* (x .+ 5) .- 3.0 * t)
# define the function for the initial values for the exact solution of the
# NLS-equation
init2(x) = exact_sul(x, 0)

function explicit_NLS2(init_val::Function, dx::Real, dt::Real; lambda::Real = -2.0,T::Real=6.0)

    # Only save the proportion-th number of iterations, 
    # this allows smaller time steps without running out of memory,
    # so only 1/proportion is of the time steps are saved
    proportion = 100

    # compute the number of dimensions for the matrices
    J = Int(30.0 / dx)
    K = Int(T / dt) * proportion
    dt = dt/proportion

    # define the size of the matrix
    m_size = 2 * J - 1

    # define matrix for the 2nd derivative
    A = spdiagm(0 => 2 * ones(m_size),
                -1 => -1 * ones(m_size - 1),
                1 => -1 * ones(m_size - 1))


    # define solution matrix
    U = zeros(ComplexF64, m_size,Int(K / proportion))

    # discretisation of the spatial variable
    xs = collect(range(start=-30, stop=30, length=m_size))

    # compute the initial values along the spatial domain
    U[:, 1] = init_val.(xs)
    Uk = deepcopy(U[:, 1])

    # Iteration loop for the time stepping
    for k in 1:K-1
        if k/proportion >= size(U, 2)
            break
        end
        # compute the nonlinear term
        Fk = 0.25 * lambda * (2 * abs2.(Uk)) .*(2 * Uk)
        # compute the next vector using the FD-scheme
        Uk1 = Uk + (dt / (abs2(dx) * im)) * A * Uk - (dt / im) * Fk

        if k < 10 || maximum(abs.(Uk1))>3
            @info "$(maximum(abs.(Uk1)))"
        end

        if k % proportion == 0
            U[:, floor(Int, k/proportion+0.01)+1] = Uk1
        end
        if k % proportion == 0 & k != 1
            Uk = deepcopy(U[:, floor(Int, k/proportion+0.01)+1])
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

dt_expl = 0.00001    # time step size for explicit scheme
dx_expl = 0.1     # mesh size for explicit scheme

@time U_num_exp = explicit_NLS2(init2, dx_expl, dt_expl, T = 0.45)
@info "Computation with the explicit scheme - done"
@info "Mesh size explicit scheme: $(Printf.@sprintf("%.4e", dx_expl))"
@info "time step explicit scheme: $(Printf.@sprintf("%.4e", dt_expl))"

# comparision of the first conservation law
# start with the explicit scheme
con1_val, con1_max, con1_min, con1_diss = first_conservation(U_num_exp)
con1_val = con1_val[.!isnan.(con1_val)] # remove any NAN values from the instability
con1_max = maximum(con1_val)
con1_min = minimum(con1_val)
con1_diss = con1_max - con1_min
@info "+++ 1st conservation law +++"
@info "The explicit scheme"
@info "Maximal value of the energy: $(Printf.@sprintf("%.4e", con1_max))"
@info "Minimal value of the energy: $(Printf.@sprintf("%.4e", con1_min))"
@info "Energy dissipation (max - min): $(Printf.@sprintf("%.4e", con1_diss))"