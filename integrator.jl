module integrator
    using input
    using TB_types
    using parms
    using oprs
    using MuladdMacro
    using ParallelAccelerator
    export integrator_t, integrator_init

    struct int_cache_t
        k1 :: TB_data_t
        k2 :: TB_data_t
        k3 :: TB_data_t
        k4 :: TB_data_t
        tmp :: TB_data_t
        dw :: TB_data_t
        u_half :: TB_data_t
    end

    struct integrator_t
        integrate :: Function
        cache :: int_cache_t
    end

    function int_cache_init(N)
        tmplt = zero_TB_data(N)
        k1 = deepcopy(tmplt)
        k2 = deepcopy(tmplt)
        k3 = deepcopy(tmplt)
        k4 = deepcopy(tmplt)
        tmp = deepcopy(tmplt)
        dw = deepcopy(tmplt)
        u_half = deepcopy(tmplt)

        return int_cache_t(k1, k2, k3, k4, tmp, dw, u_half)
    end

    @inline function runge_kutta_2nd_order(w :: TB_data_t, opr :: oprs_t, int_cache, dt)

        @muladd begin
            opr.f(w, int_cache.tmp)

            @. int_cache.u_half = w + dt * int_cache.tmp
            @. int_cache.u_half .*= opr.cache.con.ifactor

            opr.f(int_cache.u_half, int_cache.k1)

            @. int_cache.k2 = int_cache.tmp * dt/2 + w
            @. int_cache.k2 .*= opr.cache.con.ifactor

            @. w = dt/2 * int_cache.k1 + int_cache.k2
        end
        nothing
    end

    @inline function runge_kutta(w :: TB_data_t, opr :: oprs_t, int_cache, dt)

        dt_2 = 0.5 * dt

        @muladd begin
            # 1
            opr.f(w, int_cache.k1);

            # 2
            @. int_cache.tmp = w + dt_2 * int_cache.k1
            opr.f(int_cache.tmp, int_cache.k2);

            # 3
            @. int_cache.tmp = w + dt_2 * int_cache.k2
            opr.f(int_cache.tmp, int_cache.k3);

            # 4
            @. int_cache.tmp = w + int_cache.k3 * dt
            opr.f(int_cache.tmp, int_cache.k4);

            # final
            @. int_cache.dw = Complex(1/6 * dt) * (int_cache.k1+ Complex(2.0) *int_cache.k2 + Complex(2.0) * int_cache.k3 + int_cache.k4);

            @. w .+= int_cache.dw
        end

        nothing
    end

    function marching(step :: Function, pms :: parms_t, opr :: oprs_t, sol_s :: sol_s_t,
            int_cache :: int_cache_t, output_plan)

        n_last = size(sol_s.time,1)
        w = deepcopy(sol_s.sol[n_last])
        t = sol_s.time[n_last]
        max_tstep_cnt = output_plan.max_tstep_cnt

        for i = 1 : max_tstep_cnt
            @time step(w, opr, int_cache, pms.dt);

            t = t+ pms.dt
            e  = real(sum( opr.cache.tmp.μ .* conj(opr.cache.tmp.μ) )) ./ pms.N^3
            #println("max velocity", maximum(real(opr.cache.tmp.μ)))
            if mod(i, output_plan.output_freq) == 0
                    push!(sol_s.sol, deepcopy(w))
                    push!(sol_s.time, deepcopy(t))
                    push!(sol_s.k_energy, copy(e))
            end

            println("t= ", real(t), ", i= ", i)
            println("e= ", e)
            end

        return sol_s;
    end

    function integrator_init(parms)
        N = parms.N
        if parms.RK == 2
            println("Using 2 stages Runge Kutta")
            integrate = (pms, opr, sol_s, int_cache, output_plan) -> marching(runge_kutta_2nd_order, pms, opr, sol_s, int_cache, output_plan)
        else
            println("Using 4 stages Runge Kutta")
            integrate = (pms, opr, sol_s, int_cache, output_plan) -> marching(runge_kutta, pms, opr, sol_s, int_cache, output_plan)
        end

        return integrator_t(integrate, int_cache_init(N))
    end
end
