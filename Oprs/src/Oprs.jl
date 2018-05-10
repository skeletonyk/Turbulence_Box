module Oprs
    using Parms
    using TB_types
    using MuladdMacro
    export oprs_t, oprs_init

    include("opr_types.jl")
    include("curl.jl")
    include("nonlin.jl")

    function oprs_init(parms)
        cache = opr_cache_init(parms.N)
        f = fw
        return oprs_t(f, cache, parms.ν );
    end

    @inline function dot_product(a :: TB_data_t, b :: TB_data_t, c :: TB_data_t)
    for i = 1 : 3
        @muladd  @. c[i] = a[i] * b[i]
    end
    end


    @inline function fw(w :: TB_data_t, a :: TB_data_t, cache, ν)

        w[2][1,1,1].=0
        w[1][1,1,1].=0
        w[3][1,1,1].=0

        nonlin(w, cache)

        dot_product(cache.con.l_, w, cache.tmp.lap_w)
        for i = 1:3
            @muladd @. a[i] = - cache.tmp.nl[i] - ν * cache.tmp.lap_w[i]
        end
        nothing
    end

end
