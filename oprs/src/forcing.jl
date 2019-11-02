function forcing_calc(u :: TB_data_t, force :: TB_data_t, f_mag :: Float64)
    N = size(u,3)
    shell_ind   = [1,2,3,4, N, N-1,N-2]
    shell_k     = [0,1,2,3,1,2,3]

    cnt = 0
    δ = Complex(f_mag)#0.0010662224073302792#0.017059558517284468
    for kk = 1:7
        for jj = 1:7
            for ii =  1:4
            for m = 1:3
                i = shell_ind[ii]; j = shell_ind[jj]; k = shell_ind[kk];
                k_norm = sqrt(shell_k[ii]^2 + shell_k[jj]^2 +shell_k[kk]^2)
                if (k_norm<=3.5 && k_norm>=0.1)
                    #println(i," ",j," ",k)
                    if abs(u[i,j,k,m])>0
                        cnt +=1
                        force[i,j,k, m] = δ /conj(u[i,j,k,m])#*(u[i,j,k,:]) #/u[i,j,k,:] #conj(u[i,j,k,:]) / (u[i,j,k,:] * conj(u[i,j,k,:]))
                    end
                end
            end
            end
        end
    end
    #println(cnt)
    #force .*= (1/cnt/2)*N^6
    force .*= (1/cnt/2)*N^6 /(2*pi)^3

    nothing

end
