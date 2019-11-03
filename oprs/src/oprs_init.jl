using FFTW
    function oprs_init(parms, minimum=false)
        cache = opr_cache_init(parms, minimum);

        if parms.RK == 4
            f = (w,a) -> f_wo_if_exec(w,a,cache,parms.ν, parms.forcing, parms.f_mag)
        elseif parms.RK == 2
            f = (w,a) -> f_if_exec(w,a,cache,parms.ν, parms.forcing, parms.f_mag)
        else
            error("RK need to be either 2 or 4")
            nothing
        end

        curl = (x, y) -> curl_exec(x, y, cache.con.k)
        return oprs_t(f, curl, cache, parms.ν);
    end

    function opr_cache_init(parms, minimum=false)

        if (minimum)
            N = 1
            N_pad =2
            println("Not building cache temp")
        else
            N = parms.N
            N_pad = parms.N_pad
        end

        tmplt = zero_TB_data(N)
        tmplt_big = zero_TB_data(N_pad)
        tmplt_big_real = zero_TB_real_data(N_pad)
        tmplt_big_half = zeros(Complex{Float64}, N_pad>>1+1, N_pad, N_pad,3)

        tmp = opr_tmp_cache_t(
            deepcopy(tmplt),
            deepcopy(tmplt),
            deepcopy(tmplt),
            deepcopy(tmplt),
            deepcopy(tmplt),
            deepcopy(tmplt_big),
            deepcopy(tmplt_big_real),
            deepcopy(tmplt_big_real),
            deepcopy(tmplt_big_real),
            deepcopy(tmplt),
            deepcopy(tmplt)
            );

        con = opr_const_cache_build(parms, minimum)

        return opr_cache_t(tmp, con)
    end

    function opr_const_cache_build(parms, minimum=false)
        N = parms.N
        N_pad = parms.N_pad
        ν = parms.ν

        if minimum
            println("only builduing l_inv and curl")
            k      = zero_TB_data(N)
            l_     = zero_TB_data(1)
            l_inv  = zero_TB_data(N)
            lc_inv = zero_TB_data(1)
            ifactor = zero_TB_data(1)
        else
            k      = zero_TB_data(N)
            l_     = zero_TB_data(N)
            l_inv  = zero_TB_data(N)
            lc_inv = zero_TB_data(N)
            ifactor = zero_TB_data(N)
        end

        k0  = [collect(0: (N/2 -1)) ;0 ; collect(-N/2+1: -1) ] .* im;

        # for l_inv
        for z = 1:3
            for m = 1:N
                for j = 1:N
                    for i = 1:N>>1+1
                        tmp = (k0[i].^2 + k0[j].^2 +k0[m].^2)
                        if abs(tmp)>0.5
                            l_inv[i,j,m,z] = 1/(k0[i].^2 + k0[j].^2 +k0[m].^2)
                        else
                            l_inv[i,j,m,z] = 0
                        end

                    end
                end
            end
        end

        # for k
        for i = 1:N>>1+1
            k[i,:,:,1] .= k0[i]
        end
        for i = 1:N
            k[:,i,:,2] .= k0[i]
        end
        for i = 1:N
            k[:,:,i,3] .= k0[i]
        end

        if minimum
            N_2 = 1
            N_pad=2

            ind = [1:N_2 ; N_pad - N_2 + 1 : N_pad ]

            f_pad_half = zero_TB_data(N_pad)

            f_inv   = plan_irfft(f_pad_half, N_pad, [1,2,3], flags = FFTW.MEASURE)
            f_      = plan_rfft(rand(N_pad,N_pad,N_pad,3), [1,2,3], flags = FFTW.MEASURE)

            f_view  = view(f_pad_half, 1:N_2+1, ind, ind,:)

            @. ifactor = exp.(l_ .* ν .* parms.dt)
        else

            # for l_
            for z = 1:3
                for m = 1:N
                    for j = 1:N
                        for i = 1:N>>1+1
                            l_[i,j,m,z] = (k0[i]^2 + k0[j]^2 +k0[m]^2)
                        end
                    end
                end
            end

            N_2 = N>>1

            ind = [1:N_2 ; N_pad - N_2 + 1 : N_pad ]

            f_pad_half = zero_TB_data(N_pad)

            f_inv   = plan_irfft(f_pad_half, N_pad, [1,2,3], flags = FFTW.MEASURE)
            f_      = plan_rfft(rand(N_pad,N_pad,N_pad,3), [1,2,3], flags = FFTW.MEASURE)

            f_view  = view(f_pad_half, 1:N_2+1, ind, ind,:)

            @. ifactor = exp.(l_ .* ν .* parms.dt)
        end


        return opr_const_cache_t(k    ,
                                l_    ,
                                l_inv ,
                                ind   ,
                                f_    ,
                                f_inv ,
                                f_view,
                                f_pad_half,
                                ifactor)

    end
