module Integrator
    using TB_types
    using Parms
    using Oprs
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
    end

    struct integrator_t
        marching :: Function
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

        return int_cache_t(k1, k2, k3, k4, tmp, dw)
    end

    function integrator_init(N)
        return integrator_t(marching, int_cache_init(N))
    end

    @inline function runge_kutta(w :: TB_data_t, opr :: oprs_t, int_cache, dt)
        dt_2 = 0.5 * dt

        @muladd begin
        # 1
        opr.f(w, int_cache.k1, opr.cache, opr.ν);

        # 2
        for i = 1:3
            @. int_cache.tmp[i] = w[i] + dt_2 * int_cache.k1[i]
        end
        opr.f(int_cache.tmp, int_cache.k2, opr.cache, opr.ν);

        # 3
        for i = 1:3
            @. int_cache.tmp[i] = w[i] + dt_2 * int_cache.k2[i]
        end
        opr.f(int_cache.tmp, int_cache.k3, opr.cache, opr.ν);

        # 4
        for i = 1:3
            @. int_cache.tmp[i] = w[i] + int_cache.k3[i] * dt
        end
        opr.f(int_cache.tmp, int_cache.k4, opr.cache, opr.ν);

        # final
        for i = 1:3
            @. int_cache.dw[i] = Complex(1/6 * dt) * (int_cache.k1[i] + Complex(2.0) *int_cache.k2[i] + Complex(2.0) * int_cache.k3[i] + int_cache.k4[i]);
        end
        end

        nothing
    end

    function marching(pms :: parms_t, opr :: oprs_t, sol :: sol_s_t, int_cache :: int_cache_t, N_tot)

        n_last = size(sol.time,1)
        w = sol.sol[n_last]
        t = sol.time[n_last]

        for i = 1 : N_tot
            @time runge_kutta(w, opr, int_cache, pms.dt);
            w .+= int_cache.dw
            t = t+ pms.dt
            println("t= ", real(t), ", i= ", i)
        end

        nothing
    end
end
