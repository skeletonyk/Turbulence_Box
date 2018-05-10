module Parms

    export parms_init, abstract_parms_t
    export parms_t

    abstract type abstract_parms_t end

    struct parms_t{I, T} <: abstract_parms_t
        N :: I
        k :: Array{T,1}
        L :: T
        a :: T
        Re :: T
        ν :: T
        dx :: T
        dt :: T
        U_max :: T
    end

    function parms_init(N, Re)
        k = complex([collect(0: (N/2 -1)) ; collect(-N/2: -1) ]) * 2 *pi;
        Re = complex(Float64(Re))
        a = complex(0.0)
        L = complex(1.0)
        ν = 1.0/Re
        U_max = complex(1.0)
        dx = L ./ N
        dt = complex(2.8 * dx/(pi * U_max))/5

        return parms_t(N, k, L, a, Re, ν, dx, dt, U_max)
        nothing
    end

end
