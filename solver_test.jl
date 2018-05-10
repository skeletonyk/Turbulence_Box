push!(LOAD_PATH, pwd())
    FFTW.set_num_threads(Sys.CPU_CORES)
    #FFTW.set_num_threads(1)
import Solver;
reload("Solver")
N = 64
Re = 100.0
N_tot = 10

solver = Solver.init_solver(N, Re, N_tot, "KevinHelmholtz");

solver.solve(solver)
