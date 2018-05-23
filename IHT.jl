push!(LOAD_PATH, pwd())
FFTW.set_num_threads(Sys.CPU_CORES)
#FFTW.set_num_threads(1)
import input;reload("input")
import solver;reload("solver")
import sol_show;
using JLD
dir = "/home/ke/Dropbox/LES_stat/Turbulence_box"

N = 32
Re = 1/0.017/2
max_tstep_cnt = 1
output_freq = 100
ICC = "IHT"
fname = dir * "/data/" * ICC * "_N_" * string(N) * "_Re_" * string((round(Re))) * ".jld"
RK = 2
dt = 0.025
ICC = "restart"

#Forcing = false
Forcing = true

input_ = input.input_init(N, Re, ICC, max_tstep_cnt, output_freq, RK, dt, Forcing)

if ICC == "restart"
    s = load(fname)
    sols = s["sols"]
    solver_ = solver.init_solver(input_, sols);
else
    solver_ = solver.init_solver(input_);
end
sols = solver_.solve(solver_);
save(fname, "sols", sols)

import TB_stats; reload("TB_stats")
#using Plots
#gr(show = true)
w1 = sols.sol[ length(sols.sol)]
Ek = TB_stats.E_k(w1, solver_.opr)
Ek = Ek[2:end]
k = 2:size(w1,3)>>1
#plot(k, Ek, yscale =:log10, xscale =:log10)
#plot!(k, k.^(-5/3), yscale =:log10, xscale =:log10)
