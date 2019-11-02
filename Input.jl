module Input
    export input_t, input_init, output_plan_t

    struct input_t
        N :: Int
        Re :: Real
        IC :: String
        dt :: Real
        max_tstep_cnt :: Int
        output_freq :: Int
        RK :: Int
        forcing :: Bool
        f_mag :: Float64
    end

    struct output_plan_t
        max_tstep_cnt
        output_freq
    end

    function input_init(N, Re, IC, max_tstep_cnt, output_freq, RK, dt = 0.0, forcing = false, f_mag=0.0)
        return input_t(N, Re, IC, dt, max_tstep_cnt, output_freq, RK, forcing, f_mag)
    end
end
