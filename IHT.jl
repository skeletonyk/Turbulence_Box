push!(LOAD_PATH, pwd())
using FFTW
FFTW.set_num_threads(8)
using TB_types
import input;
using solver;
import sol_show;
using JLD2
dir = "./"

N = 256
Re = 300
Re = 400 * (N/32)^(4/3)
#Re = 1/0.017/2
max_tstep_cnt = 1
output_freq = 2000
ICC = "IHT"
fname = dir * "/data/" * ICC * "_N_" * string(N) * "_Re_" * string((round(Re))) * ".jld"
RK = 2
dt = 0.01
ICC = "restart"
f_mag = (1/Re)^3 *(N/2/1.4)^4
dt = sqrt(abs(f_mag)*Re)/2 * 0.1


Forcing = false
Forcing = true

input_ = input.input_init(N, Re, ICC, max_tstep_cnt, output_freq, RK, dt, Forcing, f_mag)

if ICC == "restart"
    #s = load(fname)
    @load fname sols
    #sols = s["sols"]
    solver_ = solver.init_solver(input_, sols);
else
    solver_ = solver.init_solver(input_);
end
sols = solver_.solve(solver_);
#save(fname, "sols", sols)
@save fname sols

#using Plots
#gr(show = true)
#w1 = sols.sol[ length(sols.sol)]
#Ek = TB_stats.E_k(w1, solver_.opr)
#Ek = Ek[2:end]
#k = 2:size(w1,3)>>1
#plot(k, Ek, yscale =:log10, xscale =:log10)
#plot!(k, k.^(-5/3), yscale =:log10, xscale =:log10)

using TB_stats
w = sols.sol[end]
opr = solver_.opr
stat=stats_build(w,solver_.opr,1/Re)
include("./figs/IHT_plot.jl")
println(stat.ηk_max)

Ek=stat.Ek
k=(0.5) : (size(Ek,1)-0.5)
l=sum(Ek./k)*pi/2/stat.u_prime^2
