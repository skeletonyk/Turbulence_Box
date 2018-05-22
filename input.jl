module input
    export input_t, input_init, output_plan_t

    struct input_t
        N :: Int
        Re :: Real
        IC :: String
        max_tstep_cnt :: Int
        output_freq :: Int
        RK :: Int
        forcing :: Bool
    end

    struct output_plan_t
        max_tstep_cnt
        output_freq
    end

    function input_init(N, Re, IC, max_tstep_cnt, output_freq, RK, forcing = false)
        return input_t(N, Re, IC, max_tstep_cnt, output_freq, RK, forcing)
    end
end
