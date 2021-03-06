module TB_types
    export TB_data_t,TB_data_x_t, dot_product, sol_s_t, sol_init, zero_TB_data

    TB_data_x_t = Array{Complex{Float64},3}
    TB_data_t = Array{TB_data_x_t,1}
    #TB_data_x_t{T} = Array{T,3}
    #TB_data_t{T} = Array{TB_data_x_t{T},1}

    abstract type abstract_sol_s_t end
    struct sol_s_t{R} <: abstract_sol_s_t
        sol :: Vector{TB_data_t}
        time :: Vector{R}
    end

    function zero_TB_data(N)
        data_x = complex(zeros(N,N,N))
        return ( Array( [ deepcopy(data_x), deepcopy(data_x), deepcopy(data_x) ]) )
    end

    function sol_init(init_condition :: TB_data_t, t0)
        return sol_s_t( (Vector([init_condition])),  Vector([t0]))
    end

end
