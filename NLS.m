function schrodinger_fd
    % Parameters
    lambda = 1;   % Nonlinear parameter
    epsilon = 1e-6; % Fixed-point iteration tolerance
    L = 30;       % Domain size [-L, L]
    T = 6;        % Time range [0, T]
    dx = 0.1;     % Spatial step size
    dt = 0.01;    % Temporal step size
    J = floor(L / dx); % Number of spatial points
    K = floor(T / dt); % Number of time steps

    x = -L:dx:L; % Spatial grid
    x = linspace(-L,L, 2*J-1)
    Nx = length(x);
    
    % Initial condition
    u0 = initial_condition(x);
    
    % Sparse matrix A for second derivative
    main_diag = 2 * ones(Nx, 1);
    off_diag = -1 * ones(Nx - 1, 1);
    A = spdiags([off_diag, main_diag, off_diag], -1:1, Nx - 1, Nx - 1);

    % Scaling for operator
    I = speye(Nx - 2);
    operator1 = (1i / dt) * I + (1 / (2 * dx^2)) * A;
    operator2 = (1i / dt) * I - (1 / (2 * dx^2)) * A;

    % Time-stepping loop
    U = u0(:); % Current solution, including boundary points
    U_next = U; % Placeholder for the next time step
    conservation1 = zeros(1, K);
    conservation2 = zeros(1, K);

    for k = 1:K
        U_interior = U(2:end-1); % Remove boundary values
        U_new = U_interior; % Initialize fixed-point iteration
        
        % Fixed-point iteration
        while true
            nonlinear_term = compute_nonlinear(U_new, U_interior, lambda);
            rhs = operator2 * U_interior - nonlinear_term;
            U_next_interior = operator1 \ rhs; % Solve sparse linear system
            
            if norm(U_next_interior - U_new, inf) < epsilon
                break;
            end
            U_new = U_next_interior;
        end

        % Update solution
        U_next(2:end-1) = U_next_interior;
        U_next(1) = 0; U_next(end) = 0; % Enforce boundary conditions
        U = U_next;

        % Conservation properties
        conservation1(k) = compute_conservation1(U, dx);
        conservation2(k) = compute_conservation2(U, dx, lambda);
    end

    % Plot results
    figure;
    plot(x, abs(U));
    title('Final Solution |u|');
    
    figure;
    plot(1:K, conservation1, 1:K, conservation2);
    legend('Mass Conservation', 'Energy Conservation');
    title('Conservation Properties');
end

function u0 = initial_condition(x)
    % Initial condition as per the given exact solution
    u0 = (3/2) * exp(1i * (2 * x)) .* sech((3/2) * (x + 5));
end

function F = compute_nonlinear(U_new, U, lambda)
    % Compute the nonlinear term F(U_new) for fixed-point iteration
    F = (lambda / 4) * ((abs(U_new).^2 + abs(U).^2) .* (U_new + U));
end

function c1 = compute_conservation1(U, dx)
    % Mass conservation: Integral of |u|^2
    c1 = dx * sum(abs(U).^2);
end

function c2 = compute_conservation2(U, dx, lambda)
    % Energy conservation: Integral of (|u'|^2 / 2 + λ|u|^4 / 4)
    du_dx = diff(U) / dx; % Numerical derivative
    c2 = dx * (0.5 * sum(abs(du_dx).^2) + (lambda / 4) * sum(abs(U).^4));
end
