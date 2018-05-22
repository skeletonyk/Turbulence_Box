push!(LOAD_PATH, pwd())
using JLD
import TB_output; reload("TB_output")
fname = "KH_N_256_Re_300.jld"

s = load("data/" * fname)
sols = s["sols"]

#vtk= TB_output.write_vtk_single(sols.sol[21], fname)
pvd = TB_output.write_vtk_series(sols, "pvd_KH_128")
