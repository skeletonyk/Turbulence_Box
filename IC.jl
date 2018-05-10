module IC
    export KevinHelmholtz, simple, init_f

    function init_f(N, Re, IC_name)
        if IC_name == "KevinHelmholtz"
            return KevinHelmholtz(N, 2000);
        else
            return simple(N);
        end
    end

    function KevinHelmholtz(N, scale )
        ω0 = complex(zeros(N, N, N)) ;
        ω3 = complex(zeros(N, N, N)) ;
        y = collect(0:N-1) ./N;
        ω_v = scale* ( ( sech.( scale .* (y-0.25) ) ).^2 .- ( sech.( scale .* (y .- 0.75) ) ).^2 );
        for i = 1 : N
            ω3[i,:,:] .= ω_v[i]
        end
        ω3 .+= scale/10 * (rand(N,N,N)-0.5);
        ω0 .+= scale/10 * (rand(N,N,N)-0.5);
        w3 = fft(ω3)
        w0 = fft(ω0)

        return Array([w0,w0,w3]);
    end

    function simple(N)
        ω = zeros(N, N, N, 3) ;
        y = collect(0:N-1) ./N;
        ω_v = 1./(1+sin.((1:N)/N*2pi)).^2 +  (cos.((1:N)/N*2pi))
        for i = 1 : N
            ω[i,:,:,3] .= ω_v[i]
        end
        #ω .+= scale/100 * (rand(N,N,N,3)-0.5);
        w = fft(ω,[1,2,3])
        return w;
    end

end
