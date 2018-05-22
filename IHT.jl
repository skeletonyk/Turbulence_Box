push!(LOAD_PATH, pwd())
FFTW.set_num_threads(Sys.CPU_CORES)
#FFTW.set_num_threads(1)
import input;reload("input")
import solver;reload("solver")
import sol_show;
using JLD

N = 128
Re = 100
max_tstep_cnt = 10
output_freq = 50
ICC = "IHT"
RK = 2
#ICC = "restart"

Forcing = false
#Forcing = true

input_ = input.input_init(N, Re, ICC, max_tstep_cnt, output_freq, RK, Forcing)

if ICC == "restart"
    solver_ = solver.init_solver(input_, sols);
else
    solver_ = solver.init_solver(input_);
end
sols = solver_.solve(solver_);
name = "./data/IHT_N_" * string(N) * "_Re_" * string(Re) * ".jld"
save(name, "sols", sols)

import TB_stats
#using Plots
#gr(show = true)
w1 = sols.sol[1]
Ek = TB_stats.E_k(w1, solver_.opr)
k = 1:size(w1,1)>>1
#plot(k, Ek, yscale =:log10, xscale =:log10)
