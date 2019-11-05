module TB_stats
    using TB_types
    using oprs
    using FFTW
    using Statistics
    export stats_build

    struct stats_t
        Etot
        k
        Ek
        ε
        u
        μ
        λ
        Reλ
        η
        ηk_max
        u_prime
        L
    end

    function stats_build(w :: TB_data_t, opr, ν)
        N = size(w,3);
        k = 0.5:( N - 1.5)

        a = -w .* opr.cache.con.l_inv

        u = similar(w)
        opr.curl(a, u)
        # normal fft
        μ = irfft(u,N,[1,2,3])
        u = fft(μ,[1,2,3])./N^3

        Etot = 0.5*sum(μ.^2) / N^3
        Ek = Ek_calc(w, u, opr)

        ε = ε_calc(w, u, opr, ν)
        η = (ν^3/ε)^(1/4)
        u_p = sqrt(mean(μ.^2))#irfft(u,N,[1,2,3]) # u_p_calc(w, u )
        λ = (15*ν*u_p^2/ε)^(0.5) #λ_calc(u_p, ε, ν)
        Reλ = u_p*λ/ν
        ηk_max = η * N/2

        u_p=sqrt(mean(μ.^2))
        #l = sum(Ek ./ k) / sum(Ek)*π/2/u_p^2

        L = integral_scale(Ek, u_p)

        return stats_t(
        Etot,
        k,
        Ek,
        ε,
        u,
        μ,
        λ,
        Reλ,
        η,
        ηk_max,
        u_p,
        L
        )
    end
    function integral_scale(Ek, u_p)
        k = 1:1:size(Ek,1)
        return sum(Ek./k)*pi/2/u_p^2

    end

    function λ_calc(u_p, ε, ν)
        N = size(u_p, 3)
        u_p_2 = u_p.^2
        u_px_2 = (diff(u_p,1)).^2
        return mean(u_p_2)/mean(u_px_2)
    end


    function Ek_calc(w :: TB_data_t, u :: TB_data_t,opr)
        N = size(w,3)
        u_2 = sum(u .* conj(u), dims=4)

        k0  = [collect(0: (N>>1 -1)) ;0 ; collect(- (N>>1) +1: -1) ]
        Ek = zeros(Float64, N>>1 - 1)
        bin = zeros(Float64, N>>1 - 1)

        # k: 1: -> 0<= < 1
        cnt = 0
        for k = 1: N
            for j = 1:N
                for i = 1:N
                    k_r = Int( round(sqrt(k0[i].^2 + k0[j].^2 + k0[k].^2 )+0.50000000000001) )
                    #println(i," ",j," ",k," ",k_r)
                    if (k_r < (N>>1) && k_r>0)
                        cnt +=1
                        Ek[k_r] += u_2[i,j,k,1]
                        bin[k_r] += 1
                    end
                end
            end
        end
        Ek *= 1/2
        return real(Ek)
    end

    function ε_calc(w :: TB_data_t, u :: TB_data_t, opr, ν)

        N = size(w,3)
        # normal fft
        u_2 = sum(u .* conj(u), dims = [4])
        k0  = [collect(0: (N/2 -1)) ;0 ; collect(-N/2+1: -1) ]

        s = 0.0
        for k = 1: N
            for j = 1:N
                for i = 1:N
                    k² = k0[i].^2 + k0[j].^2 + k0[k].^2
                    s += u_2[i,j,k,1] * k²
                end
            end
        end
        s *= ν
        return real(s)

    end
    function dissipation_sij_sum(u, opr)

        k = opr.cache.con.k
        S=0
        for i=1:3
            for j=1:3
                sij= 0.5*(u[:,:,:,i] .* k[:,:,:,j]+u[:,:,:,j] .* k[:,:,:,i]);
                Sij= ifft(sij)
                S += sum(Sij .* conj(Sij))
            end
        end
        return S*2*nu;

    end

end
