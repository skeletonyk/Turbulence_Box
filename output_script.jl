dir = "/home/kyu/Turbulence_Box"
push!(LOAD_PATH, dir)

FFTW.set_num_threads(Sys.CPU_CORES)
using JLD

import TB_output; reload("TB_output")
fname = "KH_N_196_Re_300_long.jld"

s = load(dir * "/data/" * fname)
sols = s["sols"]

pvd = TB_output.write_vtk_series(sols, "pvd_KH_196")
