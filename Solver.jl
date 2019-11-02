module Solver
    using Input
    using TB_types
    using parms
    using ic
    using oprs
    using integrator

    export init_solver

    struct solver_t
        solve :: Function
        parms :: parms_t
        sol_s :: sol_s_t
        opr :: oprs_t
        integrator :: integrator_t
        output_plan :: output_plan_t
    end


    function init_solver(input, solution = [], filename="null")
        parms = parms_init(input);
        opr = oprs_init(parms);
        if input.IC == "restart"
            #init_condition = solution.sol[end]
            sol = solution;
        elseif input.IC == "from_file"
            init_condition = init_from_file(parms, filename, opr)
            sol = sol_init(init_condition, 0.0);
        else
            init_condition = init_f(parms, input.IC, opr);
            sol = sol_init(init_condition, 0.0);
        end
        integrator = integrator_init(parms);
        output_plan = output_plan_t(input.max_tstep_cnt, input.output_freq);

        return solver_t(solve, parms, sol, opr, integrator, output_plan);
    end

    function solve(solver)
        integrator = solver.integrator;

        return integrator.integrate(solver.parms, solver.opr, solver.sol_s, integrator.cache, solver.output_plan);
    end

end
