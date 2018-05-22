module ic
    export KevinHelmholtz, simple, init_f

    function init_f(parms, IC_name, opr)
        N = parms.N
        Re = parms.Re

        if IC_name == "KevinHelmholtz"
            w = KevinHelmholtz(N, 500);
        elseif IC_name == "simple"
            w =  simple(N);
        elseif IC_name == "IHT"
            w = IHT(N, opr)
        elseif IC_name == "ABC"
            w = ABC(N, opr)
        end

        N_2 = size(w,3)>>1 + 1
        w[N_2, :, :, :] .= 0
        w[:, N_2, :, :] .= 0
        w[:, :, N_2, :] .= 0
        return w[1:N_2, :, :, : ]
    end

    function KevinHelmholtz(N, scale )
        Ascale = N/128.0*30
        ω0 = complex(zeros(N, N, N, 3)) ;
        y = collect(0:N-1) ./N;
        ω_v = Ascale*( ( sech.( scale .* (y-0.25) ) ).^2 .- ( sech.( scale .* (y .- 0.75) ) ).^2 );
        for i = 1 : N
            ω0[i,:,:,3] .= ω_v[i]
        end
        ω0 .+= Ascale/100 * (rand(N,N,N,3)-0.5);
        w0 = fft(ω0,[1,2,3])

        return w0
    end

    function IHT(N, opr)
        w = (zeros(N, N, N, 3)) ;
        y = collect(0:N-1) ./N;

        k0 = [collect(0: (N/2 -1)) ;0 ; collect(-N/2+1: -1) ]

        for k = 1:N
            for j = 1:N
                for i = 1:N
                    for m = 1:3
                    k_bar = sqrt(k0[i].^2 + k0[j].^2 + k0[k].^2)
                    if k_bar >0
                        w[i,j,k,m] .= rand()#w[i,j,k] = 3e4 * rand()*k_bar^(-2/3)
                    else
                        w[i,j,k,m] .= rand()
                    end
                end
                end
            end
        end
        ww = similar(rfft(w,[1,2,3]))
        opr.curl(rfft(w,[1,2,3]), ww)
        return ww
    end


    function taylor_green(N)
    end

    function ABC(N, opr)
        A = 1
        B = 1
        C = 0

        w = complex(zeros(N, N, N, 3) );
        μ = zeros(N, N, N, 3) ;

        y = collect(0:N-1) ./N * 2*pi;

        for k = 1:N
            for j = 1:N
                for i = 1:N
                    μ[i,j,k,1] = A * sin(y[k]) + C * cos(y[j])
                    μ[i,j,k,2] = B * sin(y[i]) + A * cos(y[k])
                    μ[i,j,k,3] = C * sin(y[j]) + B * cos(y[i])
                end
            end
        end

        u = fft(μ, [1,2,3])
        opr.curl(u,w)

        return w
    end

    function simple(N)
        ω = zeros(N, N, N, 3) ;
        y = collect(0:N-1) ./N;
        ω_v = 1./((1+sin.((1:N)/N*2pi)).^2 +  (cos.((1:N)/N*2pi))+1)
        for i = 1 : N
            ω[i,:,:,3] .= ω_v[i]
        end
        #ω .+= scale/100 * (rand(N,N,N,3)-0.5);
        w = fft(ω,[1,2,3])
        return w;
    end

end
