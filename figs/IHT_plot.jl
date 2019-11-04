dir = "."
push!(LOAD_PATH, dir)
using Statistics
import TB_output;
#import oprs;

#s = load(fname)
#sols = s["sols"]
w = sols.sol[ length(sols.sol)]
import TB_stats;
using Plots
#using PGFPlots
using LaTeXStrings
#PGFPlots.pushPGFPlotsPreamble("\\usepackage{amssymb}")
ν = 1/Re

stats = TB_stats.stats_build(w, opr, ν)
#gr(show = true)
Ek = stats.Ek[3:end]
#k = collect(stats.k[3:end])
k = collect(2.5:1:2.5+size(Ek,1)-1)
ηk = stats.η * k
y = Ek ./ (stats.ε * ν^5)^(1/4) .* (ηk.^(5/3))
y2 = Ek./(stats.ε * ν^5)^(1/4)

plot(ηk,
y2,
lw=2,
xlims = (0.03,1.5),
style=:auto,
yscale =:log10,
xscale =:log10)

##
plot(ηk,
y,
lw=2,
xlims = (0.03,1.5),
style=:auto,
yscale =:log10,
xscale =:log10,
label = L"$$ \eta E(k) / (\varepsilon \nu^5)^{1/4} $$",
xlabel = L"$ \eta  k $",
ylabel = L"$$ E(k) $$"
)

figname = dir * "/figs/" * "JHTDB_N_" * string(N) * "_Re_" * string((round(Re)))
savefig(figname*"Ek_flat.pdf")
savefig(figname*"Ek_flat.png")

# plot 2
plot(ηk,
y2,
lw=2,
#xlims = (0.03,1.5),
style=:auto,
yscale =:log10,
xscale =:log10,
label  = L"$\eta E(k) / (\varepsilon \nu^5)^{1/4}$ ",
xlabel = L"$\eta  k$ ",
ylabel = L"$\eta E(k)/ (\varepsilon \nu^5)^{1/4}$ "
#label = L"$$ \eta E(k) / (\varepsilon \nu^5)^{1/4} $$",
#xlabel = L"$ \eta  k $",
#ylabel = L"$$ \eta E(k) / (\varepsilon \nu^5)^{1/4} $$"
)

plot!(ηk,
stats.η * k.^(-5/3)*1000000,
#ηk.^(-5/3)*2,
lw=2,
style=:auto,
yscale =:log10,
xscale =:log10,
label = L"-5/3",
)


#xlabel!(L"$$ \eta k $$")
#ylabel!(L"$$ E(k) $$") #"E(k)")

figname = dir * "/figs/" * "JHTDB_N_" * string(N) * "_Re_" * string((round(Re)))
savefig(figname*"Ek.pdf")
savefig(figname*"Ek.png")

plot(sols.k_energy,style=:auto,label="kinetic enerygy", lw = 3, title = "Total kinetic energy with time")
xlabel!("T")
ylabel!("K")
savefig(figname*"k_vs_t.pdf")
savefig(figname*"k_vs_t.png")
