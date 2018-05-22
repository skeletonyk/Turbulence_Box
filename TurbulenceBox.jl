module TurbulenceBox
    using Input
    using Solver
    FFTW.set_num_threads(Sys.CPU_CORES)
    export runn

    function runn()
    N = 256
    Re = 100.0
    N_tot = 5

    solver = Solver.init_solver(N, Re, N_tot, "KevinHelmholtz");
    opr = solver.opr;
    conopr = opr.cache.con

    solver.solve(solver)
    end
end
