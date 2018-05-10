@inline function pad_ifft(u :: TB_data_t, v :: TB_data_t, f_inv, f_view, f_pad)
    for i = 1 : 3
        @muladd f_view[i] .= u[i]
    end

    for i = 1 : 3
        A_mul_B!(v[i], f_inv, f_pad[i])
    end
    nothing
end

@inline function cut(a :: TB_data_x_t, b :: TB_data_x_t, ind :: Array{Int,1})

    b .= view(a, ind, ind, ind) #a[ ind, ind, ind]
    #for i = 1 : length(ind)
    #    for j = 1 : length(ind)
    #        for k = 1: length(ind)
    #            @. b[k,j,i] .= a[ind[k], ind[j], ind[i]]
    #        end
    #    end
    #end

end

@inline function nonlin{F}(w :: TB_data_t, cache :: F)
    # 1.
    dot_product(w, cache.con.lc_inv, cache.tmp.u)

    # ω × u
    # get u, w in real space
    pad_ifft(cache.tmp.u, cache.tmp.μ, cache.con.f_inv, cache.con.f_view, cache.con.f_pad)
    #println(cache.μ[:,1,1,1])
    pad_ifft(   w, cache.tmp.ω, cache.con.f_inv, cache.con.f_view, cache.con.f_pad)

    # cross in real space (2N^3)
    cross_3d(cache.tmp.ω, cache.tmp.μ, cache.tmp.ωμ)

    # Fourier back
    for i = 1:3
        A_mul_B!(cache.tmp.wu[i], cache.con.f_, cache.tmp.ωμ[i])
        cut(cache.tmp.wu[i], cache.tmp.wu_cut[i], cache.con.ind)
        #cache.tmp.wu_cut[i] .= cache.tmp.wu[i][ cache.con.ind, cache.con.ind, cache.con.ind]
    end

    curl(cache.tmp.wu_cut, cache.tmp.nl, cache.con.k)
    nothing
end
