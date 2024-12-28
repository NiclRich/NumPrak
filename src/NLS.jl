using LinearAlgebra
using Printf
include("NLS_module.jl")
include("NLSutils.jl")

using .Main.NLS
using .Main.NLSutils

# define the parameters for solving the equation
lambda = -2.0       # equation parameter

dt_dfp = 0.01       # time step size for DFP-scheme
dx_dfp = 0.1        # mesh size for DFP scheme
T_dfp  = 6.0        # end time for DFP scheme

# the time step and mesh size need to be smaller for the
# explicit scheme. It is highly unstable
dt_expl = 0.0001    # time step size for explicit scheme
dx_expl = 0.001     # mesh size for explicit scheme

# compute the numerical solution using the DFP-scheme
@time U_num_DFP = schroedingereq(init, dx_dfp, dt_dfp, T=T_dfp)
@info "Computation with the DFP scheme - done"
# compute the explicit time step scheme
@time U_num_exp = explicit_NLS(init, dx_expl, dt_expl, T = 0.45)
@info "Computation with the explicit scheme - done"

# compute the exact solution
ts = collect(range(start=0, stop=6.0, length=Int(T_dfp / dt_dfp)+1))
xs = collect(range(start=-30, stop=30.0, length=Int(30 / dx_dfp) * 2 - 1))
xs = reshape(xs, length(xs), 1)
ts = reshape(ts, 1, length(ts))
U_true = exact_sol.(xs, ts)

# compute the error in the max-norm for the DFP scheme
error_dfp = norm(abs.(U_true - U_num_DFP), Inf)
# compute the difference of the norms
error_norms = norm(abs.(U_true) - abs.(U_num_DFP), Inf)

@info "Error in the max-norm for the DFP scheme: $(Printf.@sprintf("%.4e", error_dfp))"
@info "Difference of the norms, |U_num| - |U_exact|: $(Printf.@sprintf("%.4e", error_norms))"

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

# continue with the DFP scheme
con1_val, con1_max, con1_min, con1_diss = first_conservation(U_num_DFP)
@info "The DFP scheme"
@info "Maximal value of the energy: $(Printf.@sprintf("%.4e", con1_max))"
@info "Minimal value of the energy: $(Printf.@sprintf("%.4e", con1_min))"
@info "Energy dissipation (max - min): $(Printf.@sprintf("%.4e", con1_diss))"

# comparision of the second conservation law
# start with the explicit scheme
con1_val, con1_max, con1_min, con1_diss = second_conservation(U_num_exp, dx_expl, lambda)
con1_val = con1_val[.!isnan.(con1_val)] # remove any NAN values from the instability
con1_max = maximum(con1_val)
con1_min = minimum(con1_val)
con1_diss = con1_max - con1_min
@info "+++ 2nd conservation law +++"
@info "The explicit scheme"
@info "Maximal value of the energy: $(Printf.@sprintf("%.4e", con1_max))"
@info "Minimal value of the energy: $(Printf.@sprintf("%.4e", con1_min))"
@info "Energy dissipation (max - min): $(Printf.@sprintf("%.4e", con1_diss))"

# continue with the DFP scheme
con1_val, con1_max, con1_min, con1_diss = second_conservation(U_num_DFP, dx_dfp, lambda)
@info "The DFP scheme"
@info "Maximal value of the energy: $(Printf.@sprintf("%.4e", con1_max))"
@info "Minimal value of the energy: $(Printf.@sprintf("%.4e", con1_min))"
@info "Energy dissipation (max - min): $(Printf.@sprintf("%.4e", con1_diss))"

# create gifs 
create_gif(U_num_DFP, U_true, "Absolute value of solution", dt_dfp, T_dfp, "NLS_abs.gif"; projection = "abs")
create_gif(U_num_DFP, U_true, "Real part of solution", dt_dfp, T_dfp, "NLS_real.gif"; projection = "real")
create_gif(U_num_DFP, U_true, "Imaginary part of solution", dt_dfp, T_dfp, "NLS_imag.gif"; projection = "imag")