function forcing_calc(u :: TB_data_t, force :: TB_data_t, f_mag :: Float64)
    N = size(u,3)
    shell_ind   = [1,2,3,4, N, N-1,N-2]
    shell_k     = [0,1,2,3,1,2,3]

    cnt = 0
    δ = Complex(f_mag)#0.0010662224073302792#0.017059558517284468

    sum_input = 0
    @. force = force *0;
    for kk = 1:7
        for jj = 1:7
            for ii = 1:4
            for m = 1:3
                i = shell_ind[ii];
                j = shell_ind[jj];
                k = shell_ind[kk];
                k_norm = sqrt(shell_k[ii]^2 + shell_k[jj]^2 +shell_k[kk]^2)

                if (k_norm<2.5 && k_norm>=0.01)
                    #println(i," ",j," ",k)
                    #println(u[i,j,k,m])
                    if abs(u[i,j,k,m])>0.1
                        cnt +=1
                        force[i,j,k, m] = k_norm^(-8/3) / conj(u[i,j,k,m])
                        if (shell_k[ii] == 0)
                            sum_input += k_norm^(-8/3)
                        else
                            sum_input += k_norm^(-8/3)*2
                        end
                    end
                end
            end
            end
        end
    end
    force .*= (f_mag/sum_input)*N^6 /(2*pi)^3

    #force .*= (1/cnt/2)*N^6 /(2*pi)^3

    #test
    #F = irfft(force, N, [1,2,3]);
    #mu =  irfft(u, N, [1,2,3]);
    #println("------------------")
    #println(sum(F .* mu) /N^3 );
    #println(sum( fft(F,[1,2,3]) .* conj(fft(mu,[1,2,3])) ) /N^6 );
    #println("------------------")

    nothing

end
