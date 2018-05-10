module Solver
    using TB_types
    using Parms
    using IC
    using Oprs
    using Integrator

    export init_solver

    mutable struct solver_t
        solve :: Function
        parms :: parms_t
        sol_s :: sol_s_t
        opr :: oprs_t
        integrator :: integrator_t
        N_tot :: Int
    end

    function init_solver(N, Re, N_tot, IC_name)
        parms = parms_init(N, Re)
        init_condition = init_f(parms.N, parms.Re, IC_name);
        sol = sol_init(init_condition, 0.0);
        opr = oprs_init(parms);
        integrator = integrator_init(parms.N)

        return solver_t(solve, parms, sol, opr, integrator, N_tot);
    end

    function solve(solver)
        integrator = solver.integrator
        integrator.marching(solver.parms, solver.opr, solver.sol_s, integrator.cache, solver.N_tot)
    end

end
