push!(LOAD_PATH, pwd())
using Revise
using Suppressor
using FFTW
FFTW.set_num_threads(16)
using TB_types
import Input;
using Solver;
import sol_show;
using JLD2
dir = "./"
f_mag = 0.103 * (2*pi)^3 #(1/Re)^3 *(N/2/1.4)^4

N = 64
nu = ( (1.5/N*2.0)*(f_mag/(2*pi)^3)^0.25 ) ^(4/3)
Re = 1/nu

Forcing = false
Forcing = true

max_tstep_cnt = 1001
output_freq = 1000
RK =    2
dt = 1/N*0.2*2*pi

ICC = "IHT"
#ICC = "from_file"
filename  = "../JHTBD/data_64.h5"

fname = dir * "data/" * ICC * "_N_" * string(N) * "_Re_" * string((round(Re))) * ".jld"

# --- generate IC from julia results of a different resolution
N_ref = 32
nu_ref = ( (1.5/N_ref*2.0)*(f_mag/(2*pi)^3)^0.25 ) ^(4/3)
Re_ref = 1/nu_ref
fname_ref = dir * "data/" * ICC * "_N_" * string(N_ref) * "_Re_" * string((round(Re_ref))) * ".jld"
# --- generate IC from julia results of a different resolution

ICC = "restart"
from_julia_file = true

#from_julia_file = false

input_ = Input.input_init(N, Re, ICC, max_tstep_cnt, output_freq, RK, dt, Forcing, f_mag);
continue_sols=true
begin
    if ICC == "restart"
        if (!from_julia_file)
            @load fname sols;
            global Solver_ = Solver.init_solver(input_, sols, "", continue_sols );
        else
            continue_sols=false
            @load fname_ref sols;
            Solver_ = Solver.init_solver(input_, sols, "", continue_sols )
        end
    elseif ICC == "from_file"
        Solver_ = Solver.init_solver(input_, [], filename)
    else
        Solver_ = Solver.init_solver(input_);
    end
end
sols = Solver_.solve(Solver_);
#save(fname, "sols", sols)
@save fname sols

#using Plots
#gr(show = true)
#w1 = sols.sol[ length(sols.sol)]
#Ek = TB_stats.E_k(w1, Solver_.opr)
#Ek = Ek[2:end]
#k = 2:size(w1,3)>>1
#plot(k, Ek, yscale =:log10, xscale =:log10)
#plot!(k, k.^(-5/3), yscale =:log10, xscale =:log10)

using TB_stats
w = sols.sol[end]
opr = Solver_.opr
stat=stats_build(w, Solver_.opr,1/Re)
#include("./figs/IHT_plot.jl")
println(stat.ηk_max)

Ek=stat.Ek
k=(0.5) : (size(Ek,1)-0.5)
l=sum(Ek./k)*pi/2/stat.u_prime^2

function dissipation_sij_sum(w, opr, nu)
        N = size(w,3);
        a = -w .* opr.cache.con.l_inv

        u = similar(w)
        opr.curl(a, u)
        # normal fft
        #μ = ifft(u,N,[1,2,3])
        #u = fft(μ,[1,2,3])./N^3
        k = opr.cache.con.k
        S=0
        for i=1:3
            for j=1:3
                sij= 0.5*(u[:,:,:,i] .* k[:,:,:,j]+u[:,:,:,j] .* k[:,:,:,i]);
                Sij= ifft(sij)
                S += sum(Sij .* conj(Sij))
            end
        end
        return S*2*nu/N^3;

end
s=dissipation_sij_sum(w,opr,nu)
print(s)
print(stat.Reλ)
