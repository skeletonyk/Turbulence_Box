push!(LOAD_PATH, pwd())
FFTW.set_num_threads(Sys.CPU_CORES)
#FFTW.set_num_threads(1)
import input;reload("input")
import solver;reload("solver")
import sol_show;
using JLD

N = 196
Re = 200
max_tstep_cnt = 1800
output_freq = 30
ICC = "KevinHelmholtz"
#ICC = "ABC"
RK = 2
#ICC = "restart"


input_ = input.input_init(N, Re, ICC, max_tstep_cnt, output_freq, RK)

if ICC == "restart"
    solver_ = solver.init_solver(input_, sols);
else
    solver_ = solver.init_solver(input_);
end
sols = solver_.solve(solver_);

fname = "./data/KH_N_" * string(N) *"_Re_" * string(Re) * ".jld"
save(fname, "sols", sols)
