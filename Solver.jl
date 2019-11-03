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


    function init_solver(input, solution = [], filename="null", continue_sols=true)
        parms = parms_init(input);
        opr = oprs_init(parms);

        if input.IC == "restart"
            #init_condition = solution.sol[end]
            if (continue_sols)
                sol = solution;
            else
                ic_ref = solution.sol[end]
                N_ref = size(ic_ref,3)
                N = input.N

                N_shared = min(N, N_ref)
                init_condition = zeros(Complex{Float64}, N>>1+1, N, N, 3)
                k_map = [collect(1:N_ref >>1+1); collect(zeros(Int64, N-N_ref))  ;N_ref-(N_ref >>1)+2:N_ref]


                for i=1: N>>1+1
                    for j = 1: N
                        for k = 1:N
                            if (k_map[i] != 0  && k_map[j] != 0  && k_map[k] != 0)
                                init_condition[i,j,k,:] .= ic_ref[k_map[i],k_map[j],k_map[k],:] * (N^3 / N_ref^3)
                            end
                        end
                    end
                end
                sol = sol_init(init_condition, 0.0);
            end
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
