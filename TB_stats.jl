module TB_stats
    using TB_types
    using oprs

    export E_k

    function E_k(w :: TB_data_t, opr)
        N = size(w,3)
        u = similar(w)
        opr.curl(w, u)
        u .*= -opr.cache.con.l_inv

        # normal fft
        u = fft(irfft(u,N,[1,2,3]),[1,2,3])./N^3

        u_2 = sum(u .* conj(u), 4)

        k0  = [collect(0: (N/2 -1)) ;0 ; collect(-N/2+1: -1) ]
        Ek = zeros(Float64, N>>1)

        for k = 1: N
            for j = 1:N
                for i = 1:N
                    k_r = Int(round(sqrt(k0[i].^2 + k0[j].^2 + k0[k].^2)))
                    if k_r < (N>>1)
                        Ek[k_r+1] .+= u_2[i,j,k,1]
                    end
                end
            end
        end

        return Ek


    end
end
