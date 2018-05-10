    struct opr_tmp_cache_t
        u ::  TB_data_t
        nl ::  TB_data_t
        wu_cut ::  TB_data_t
        lap_w ::  TB_data_t
        wu ::  TB_data_t
        μ ::  TB_data_t
        ω ::  TB_data_t
        ωμ ::  TB_data_t
    end

    struct opr_const_cache_t
        k       :: TB_data_t
        l_      :: TB_data_t
        l_inv   :: TB_data_t
        lc_inv  :: TB_data_t
        ind     :: Vector{Int64}
        f_      :: Base.DFT.FFTW.cFFTWPlan{Complex{Float64},-1,false,3}
        f_inv   :: Base.DFT.ScaledPlan{Complex{Float64},Base.DFT.FFTW.cFFTWPlan{Complex{Float64},1,false,3},Float64}
        f_view  :: Array{SubArray{Complex{Float64},3,Array{Complex{Float64},3},Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}},false},1}
        f_pad   :: TB_data_t
    end

    struct opr_cache_t
        tmp :: opr_tmp_cache_t
        con :: opr_const_cache_t
    end

    struct oprs_t
        f :: Function
        cache :: opr_cache_t
        ν
    end

    function opr_cache_init(N)
        tmplt = zero_TB_data(N)
        tmplt_big = zero_TB_data(Int(N/2*3))

        tmp = opr_tmp_cache_t(
            deepcopy(tmplt),
            deepcopy(tmplt),
            deepcopy(tmplt),
            deepcopy(tmplt),
            deepcopy(tmplt_big),
            deepcopy(tmplt_big),
            deepcopy(tmplt_big),
            deepcopy(tmplt_big)
            );

        const con = opr_const_cache_build(N)

        return opr_cache_t(tmp, con)
    end

    function opr_const_cache_build(N)
        k      = zero_TB_data(N)
        l_     = zero_TB_data(N)
        l_inv  = zero_TB_data(N)
        lc_inv = zero_TB_data(N)

        k0  = [collect(0: (N/2 -1)) ; collect(-N/2: -1) ];


        # for k
        for i = 1:N
            k[1][i,:,:] .= k0[i]
        end
        for i = 1:N
            k[2][:,i,:] .= k0[i]
        end
        for i = 1:N
            k[3][:,:,i] .= k0[i]
        end

        # for l_
        for z = 1:3
            for m = 1:N
                for j = 1:N
                    for i = 1:N
                        l_[z][i,j,m] .= (k0[i].^2 + k0[j].^2 +k0[m].^2)
                    end
                end
            end
        end

        # for l_inv
        for z = 1:3
            for m = 1:N
                for j = 1:N
                    for i = 1:N
                        l_inv[z][i,j,m] .= 1./(k0[i].^2 + k0[j].^2 +k0[m].^2)
                    end
                end
            end
            l_inv[z][1,1,1] .= 0
        end

        lc_inv = [ l_inv[i] .* k[i] for i= 1:3 ]


        N_2 = Int(N/2)

        ind = [1:N_2 ; (N+1):(N_2*3) ]
        f_pad   = zero_TB_data(N_2*3)

        f_      = plan_fft(f_pad[1], flags = FFTW.MEASURE)
        f_inv   = plan_ifft(f_pad[1], flags = FFTW.MEASURE)

        f_view  = [ view(f_pad[i],  [(1 : N_2); (N + 1) : (N_2 * 3)],
                                        [(1 : N_2); (N + 1) : (N_2 * 3)],
                                        [(1 : N_2); (N + 1) : (N_2 * 3)]) for i = 1:3]
        return opr_const_cache_t(k    ,
                                l_    ,
                                l_inv ,
                                lc_inv,
                                ind   ,
                                f_    ,
                                f_inv ,
                                f_view,
                                f_pad)

    end
