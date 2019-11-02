    struct opr_tmp_cache_t
        a ::  TB_data_t
        u ::  TB_data_t
        nl ::  TB_data_t
        wu_cut ::  TB_data_t
        lap_w ::  TB_data_t
        wu ::  TB_data_t
        μ ::  TB_data_real_t
        ω ::  TB_data_real_t
        ωμ ::  TB_data_real_t
        force :: TB_data_t
        curl_force :: TB_data_t
    end

    struct opr_const_cache_t
        k       :: TB_data_t
        l_      :: TB_data_t
        l_inv   :: TB_data_t
        ind     :: Vector{Int64}
        f_
        f_inv
        f_view
        f_pad_half :: TB_data_t
        ifactor :: TB_data_t
    end

    struct opr_cache_t
        tmp :: opr_tmp_cache_t
        con :: opr_const_cache_t
    end

    struct oprs_t
        f :: Function
        curl :: Function
        cache :: opr_cache_t
        ν
    end
