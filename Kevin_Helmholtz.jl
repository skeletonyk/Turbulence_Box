push!(LOAD_PATH, pwd())
dir = "/home/kyu/Turbulence_Box"
push!(LOAD_PATH, dir)
FFTW.set_num_threads(Sys.CPU_CORES)
import input;reload("input")
import solver;reload("solver")
import TB_types
import sol_show;
using JLD

N = 196
Re = 300
max_tstep_cnt = 3600
output_freq = 20
ICC = "KevinHelmholtz"
#ICC = "ABC"
RK = 2
fname = dir*"/data/KH_N_" * string(N) *"_Re_" * string(Re) * "_long.jld"

#ICC = "restart"

input_ = input.input_init(N, Re, ICC, max_tstep_cnt, output_freq, RK)

if ICC == "restart"
    s = load(fname)
    sols = s["sols"]
    solver_ = solver.init_solver(input_, sols);
else
    solver_ = solver.init_solver(input_);
end
sols = solver_.solve(solver_);

save(fname, "sols", sols)
