module NLSutils
using Plots
using LinearAlgebra


export exact_sol, init, create_gif
# define the explicit solution to the NLS equation
exact_sol(x,t) = 1.5 .* exp.(im .* (1.25 .* t - x )) .* sech.(1.5 .* (x .+ 5) .- 3.0 * t)
# define the function for the initial values for the exact solution of the
# NLS-equation
init(x) = exact_sol(x, 0)

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

    gif(anim, filename)
end

# create_gif(rand(10, 10), rand(10, 10), "Test", 0.1, 1.0)

end