push!(LOAD_PATH, pwd())
using JLD
import TB_output; reload("TB_output")
dir = "/home/kyu/Turbulence_Box"
fname = "KH_N_256_Re_300.jld"

s = load(dir * "/data/" * fname)
sols = s["sols"]

pvd = TB_output.write_vtk_series(sols, "pvd_KH_128")
