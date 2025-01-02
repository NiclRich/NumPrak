module NLSutils
using Plots
using LinearAlgebra


export exact_sol, init, create_gif
# define the explicit solution to the NLS equation
exact_sol(x,t) = 1.5 .* exp.(im .* (1.25 .* t - x )) .* sech.(1.5 .* (x .+ 5) .- 3.0 * t)
# define the function for the initial values for the exact solution of the
# NLS-equation
init(x) = exact_sol(x, 0)

"""
    create_gif(U_num, U_true, title::String, dt::Real, T::Real, filename::String; projection = "abs")

Creates an animated GIF that compares a numerical solution to an exact solution over time.

# Arguments
- `U_num`: A matrix representing the numerical solution. Each column corresponds to a snapshot in time, and each row represents a spatial point.
- `U_true`: A matrix representing the exact solution, structured similarly to `U_num`.
- `title::String`: The title of the animation.
- `dt::Real`: The time step size used in the simulation.
- `T::Real`: The total simulation time.
- `filename::String`: The name of the output GIF file (e.g., `"output.gif"`).

# Keyword Arguments
- `projection::String`: Specifies how to project the data for visualization. Accepted values are:
  - `"abs"`: Absolute value of the solution (default).
  - `"real"`: Real part of the solution.
  - `"imag"`: Imaginary part of the solution.

# Behavior
- The x-axis corresponds to spatial points ranging from -30 to 30.
- The y-axis corresponds to the solution values (`U_num` and `U_true`).
- At each time step, the plot compares the numerical solution (`U_num`) with the exact solution (`U_true`).
- A dynamic annotation shows the current time step (`t`) on the plot.
- The y-axis limits are adjusted based on the specified `projection`.

# Output
- Generates and saves an animated GIF to the specified `filename`.

# Example
```julia
U_num = randn(100, 50)  # Simulated numerical solution (100 spatial points, 50 time steps)
U_true = randn(100, 50)  # Simulated exact solution
title = "Numerical vs. Exact Solution"
dt = 0.1
T = 5.0
filename = "solution_comparison.gif"

create_gif(U_num, U_true, title, dt, T, filename, projection="real")
```

# Notes
Ensure the `Plots.jl` package is installed and `@animate` macro is available.
Valid projections are enforced; an error is thrown if an unknown value is provided.
"""
function create_gif(U_num, U_true, title::String, dt::Real, T::Real, filename::String; projection = "abs")
    x_label = "x"
    y_label = "u(x,t)"

    x_axis = collect(range(start=-30, stop=30, length = size(U_num, 1)))
    t_axis = collect(range(start=0.0, stop=T, length = size(U_num, 2)))

    if projection == "abs"
        y_max = max(maximum(abs.(U_num)), maximum(abs.(U_true))) + 0.5
        y_min = -0.1
    elseif projection == "real"
        y_max = max(maximum(real.(U_num)), maximum(real.(U_true))) + 0.5
        y_min = min(minimum(real.(U_num)), minimum(real.(U_true))) - 0.1
    elseif projection == "imag"
        y_max = max(maximum(imag.(U_num)), maximum(imag.(U_true))) + 0.5
        y_min = min(minimum(imag.(U_num)), minimum(imag.(U_true))) - 0.1
    else
        @error "Unknown parameter value for projection!"
    end

    
   anim =  @animate for i in 1:size(U_num, 2)
        # Extract columns for the current frame
        if projection == "abs"
            numerical_data = abs.(U_num[:, i])
            exact_data = abs.(U_true[:, i])
        elseif projection == "real"
            numerical_data = real.(U_num[:, i])
            exact_data = real.(U_true[:, i])
        elseif projection == "imag"
            numerical_data = imag.(U_num[:, i])
            exact_data = imag.(U_true[:, i])
        end
        
        # Compute a numerical value to display (example: sum of x_data)
        time = t_axis[i]
        
        # Create the plot
        p = plot(x_axis, numerical_data, 
                 xlabel=x_label, 
                 ylabel=y_label, 
                 title=title,
                 label = "Numerical solution",
                 legend=:outerright,
                 xlimits=(-30.0, 30.0),
                 ylimits=(y_min, y_max),
                 grid=true)

        plot!(x_axis, exact_data,
        label = "Exact solution")
        
        # Add the numerical value as an annotation in the top-right corner
        annotate!(maximum(x_axis), y_max - 0.3, 
                  text("t = $(round(time, digits=2))", :right, 10))
    end

    gif(anim, filename);
end
end