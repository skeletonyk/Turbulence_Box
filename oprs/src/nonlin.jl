@inline function pad_ifft(u :: TB_data_t, v :: TB_data_real_t, f_inv, f_view, f_pad_half)

    @muladd f_view .= u
    A_mul_B!(v, f_inv, f_pad_half)
    # A_mul_B! with irfft would change the value of f_pad_half somehow
    fill!(f_pad_half, complex(0.0))
    #f_pad_half .= complex(0.0)
    v .*= 1.5^3
    nothing

end

@inline function pad_fft(ωμ :: TB_data_real_t, wu_cut :: TB_data_t, cache)
    A_mul_B!(cache.tmp.wu, cache.con.f_, ωμ)

    cut(cache.tmp.wu, wu_cut, cache.con.ind)
    nothing
end

@inline function cut(a :: TB_data_t, b :: TB_data_t, ind :: Array{Int,1})

    N = size(b,3) >> 1 + 1

    b .= view(a, 1:N, ind, ind, :) #a[ ind, ind, ind]
    b[N,:,:,:] .=0
    b[:,N,:,:] .=0
    b[:,:,N,:] .=0

end

@inline function nonlin{F}(w :: TB_data_t, cache :: F)
    # 1.
    # @. cache.tmp.u = cache.con.lc_inv .* w / wrong
    curl_exec(w, cache.tmp.u, cache.con.k)
    cache.tmp.u .*= cache.con.l_inv

    # ω × u
    # get u, w in real space
    pad_ifft(cache.tmp.u, cache.tmp.μ, cache.con.f_inv, cache.con.f_view, cache.con.f_pad_half)
    pad_ifft(   w, cache.tmp.ω, cache.con.f_inv, cache.con.f_view, cache.con.f_pad_half)

    # cross in real space (2N^3)
    cross_3d(cache.tmp.ω, cache.tmp.μ, cache.tmp.ωμ)

    # Fourier back
    pad_fft(cache.tmp.ωμ, cache.tmp.wu_cut, cache)

    curl_exec(cache.tmp.wu_cut, cache.tmp.nl, cache.con.k)
    nothing
end
