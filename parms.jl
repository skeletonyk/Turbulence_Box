module parms

    export parms_init, abstract_parms_t
    export parms_t

    abstract type abstract_parms_t end

    struct parms_t{I, T} <: abstract_parms_t
        N :: I
        N_pad :: I
        k :: Array{T,1}
        L :: T
        a :: T
        Re :: T
        ν :: T
        dx :: T
        dt :: T
        U_max :: T
        RK :: Int
        forcing :: Bool
    end

    function parms_init(input)
        N = input.N
        Re = input.Re
        RK = input.RK

        N_pad = (N>>1) * 3
        k = complex([collect(0: (N/2 -1)) ; collect(-N/2: -1) ]) * 2 *pi;
        Re = complex(Float64(Re))
        a = complex(0.0)
        L = complex(1.0)
        ν = 1.0/Re
        U_max = complex(1.0)
        dx = L ./ N

        if input.dt == 0.0
            dt = complex(min(real(2.8 * dx/(pi * U_max)), real(2.8 *Re *dx^2 /pi^2)))
        else
            dt = complex(input.dt)
        end

        return parms_t(N, N_pad, k, L, a, Re, ν, dx, dt, U_max, RK, input.forcing)
        nothing
    end

end
