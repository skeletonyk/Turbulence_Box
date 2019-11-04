@inline function f_wo_if_exec(w :: TB_data_t, a :: TB_data_t, cache, ν, forcing, f_mag)
    N_2 = size(w,3)>>1 + 1

    #w[N_2, :, :, :] .= 0
    #w[:, N_2, :, :] .= 0
    #w[:, :, N_2, :] .= 0
    #w[1,1,1,:] .= 0

    nonlin(w, cache)

    @.  cache.tmp.lap_w = cache.con.l_ * w
    @muladd @. a = - cache.tmp.nl + ν * cache.tmp.lap_w

    if forcing
        cache.tmp.curl_force .=0
        forcing_calc(cache.tmp.u, cache.tmp.force, f_mag)
        curl_exec(cache.tmp.force, cache.tmp.curl_force, cache.con.k)
        @. a +=  cache.tmp.curl_force
    end

    nothing
end

@inline function f_if_exec(w :: TB_data_t, a :: TB_data_t, cache, ν, forcing, f_mag)
    N_2 = size(w,3)>>1 + 1

    w[N_2, :, :, :] .= 0
    w[:, N_2, :, :] .= 0
    w[:, :, N_2, :] .= 0

    #w[1,1,1,:] .= 0

    nonlin(w, cache)

    if forcing
        cache.tmp.curl_force .=0
        forcing_calc(cache.tmp.u, cache.tmp.force, f_mag)
        curl_exec(cache.tmp.force, cache.tmp.curl_force, cache.con.k)
        @. a = - cache.tmp.nl + cache.tmp.curl_force
    else
        @. a = - cache.tmp.nl
    end

    nothing
end
