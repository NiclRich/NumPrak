using SparseArrays

"""
    schroedingereq(dx, dt, lambda, phi, epsilon)

Sumilate the Schroedinger Equation with the constant lambda and initial value phi.

The qualtisation of time and space is given by dt and dx. epsilon is used as a termination threshold.
"""
function schroedingereq(dx, dt, lambda, phi, epsilon)
    J = Int(30/dx)
    K = Int(6/dt)

    #def matrices
    A = spdiagm(0 => 2*ones(2J-1), -1 => -1*ones(2J-2), 1 => -1*ones(2J-2))
    I = sparse(1:2J-1,1:2J-1,ones(2J-1))
    M = im/dt*I-1/(2*dx*dx)*A
    N = im/dt*I+1/(2*dx*dx)*A

    #initialise time-space array
    U = Array{ComplexF64}(undef, K+1, 2J-1)

    #initialise the initial condition
    for j in eachindex(U[1,:])
        U[1,j] = phi((j-J)*dx)
    end

    #solve the non-linear System for each time step
    for k in 2:K+1
        U[k,:] = solve_nonlineareq(U[k-1,:], lambda, epsilon, N, M)
    end

    return U
end

"""
    solve_nonlineareq(Uk :: Vector{ComplexF64}, lambda, epsilon, N, M)

Solves the equation x = M^-1 (N * y - f(x,y)) with the non-linear function f using a fixed point iteration.

M and N are matrices of appropriate dimensions for the vectors x and y of equal lenght.
The termination condition is |x_p-x_(p+1)| < epsilon, where x_p is x in itteration step p.
"""
function solve_nonlineareq(Uk :: Vector{ComplexF64}, lambda, epsilon, N, M)
    U1 = zeros(ComplexF64, size(Uk))

    while true
        #itteration setp
        U2 = M \ (N * Uk - f(Uk,U1, lambda))

        #check termination condition
        if sum(abs2,(U1 - U2)) < epsilon*epsilon
            break
        end

        U1 = U2
    end

    return U1
end

#implementation of the non-linear part from the System
function f(Uk :: Vector{ComplexF64}, U :: Vector{ComplexF64}, lambda)
    sol = similar(U)

    for j in eachindex(U)
        sol[j] = (lambda/4) * (abs2(U[j]) + abs2(Uk[j])) * (U[j] + Uk[j])
    end

    return sol
end


dx = 1              #30/dx should be an integer
dt = 0.1              #6/dt should be an integer
lambda = 1          #constant for Schrödinger Equation
epsilon = 0.0001         #termination criterion for solving of non-linear equation

#starting value for time step 0
#A - amplitude, x0 - offset, sigma - width of packet, k - frequency
function gaussian_wave_packet(x, A=1, x0=0.0, sigma=5.0, k=1.0)
    return A * exp.(-((x .- x0).^2) ./ (2 * sigma^2)) .* exp.(im * k .* x)
end

U = schroedingereq(dx, dt, lambda, gaussian_wave_packet, epsilon)

