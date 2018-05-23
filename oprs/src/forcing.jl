function forcing_calc(u :: TB_data_t, force :: TB_data_t)
    N = size(u,3)
    shell_ind   = [1,2,3,4, N, N-1,N-2]
    shell_k     = [0,1,2,3,1,2,3]

    cnt = 0
    δ = 500
    for kk = 1:7
        for jj = 1:7
            for ii =  1:4
            for m = 1:3
                i = shell_ind[ii]; j = shell_ind[jj]; k = shell_ind[kk];
                k_norm = sqrt(shell_k[ii]^2 + shell_k[jj]^2 +shell_k[kk]^2)
                if (k_norm<=2.5 && k_norm>=0.5)
                    #println(i," ",j," ",k)
                    cnt +=1
                    if abs(u[i,j,k,m])>0
                        @. force[i,j,k, m] = δ /conj(u[i,j,k,m])#*(u[i,j,k,:]) #/u[i,j,k,:] #conj(u[i,j,k,:]) / (u[i,j,k,:] * conj(u[i,j,k,:]))
                    end
                end
            end
            end
        end
    end
    #println(cnt)
    force .*= (1/cnt)*N^3

    nothing

end
