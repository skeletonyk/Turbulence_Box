push!(LOAD_PATH, pwd())
using JLD
import TB_output; reload("TB_output")
dir = "/home/ke/Dropbox/LES_stat/Turbulence_box"
fname = "KH_N_32_Re_300_long.jld"

s = load(dir * "/data/" * fname)
sols = s["sols"]

pvd = TB_output.write_vtk_series(sols, "pvd_KH_128")
