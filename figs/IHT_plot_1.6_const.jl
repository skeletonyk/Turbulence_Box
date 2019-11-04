dir = "."
push!(LOAD_PATH, dir)
using Statistics
import TB_output;

w = sols.sol[ length(sols.sol)]
import TB_stats;
using Plots
using LaTeXStrings
ν = 1/Re
stats = TB_stats.stats_build(w, opr, ν)

Ek = stats.Ek[3:end]
k = collect(3:1:3+size(Ek,1)-1)
ηk = stats.η * k

plot(k,
Ek,
lw=2,
style=:auto,
yscale =:log10,
xscale =:log10)

plot!(k,
1.6*stats.ε^(2/3)*k.^(-5/3),
lw=2,
style=:auto,
yscale =:log10,
xscale =:log10,
label = L"$ 1.6 \epsilon^{2/3} k^{-5/3}  $",
xlabel = L"$  k $",
ylabel = L"$$ E(k) $$"
)

figname = dir * "/figs/" * "JHTDB_N_" * string(N) * "_Re_" * string((round(Re))) *"const_1.6"
savefig(figname*"Ek_1.6.pdf")
savefig(figname*"Ek_1.6.png")
