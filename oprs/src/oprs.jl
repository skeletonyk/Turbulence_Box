module oprs
    using parms
    using TB_types
    using MuladdMacro
    export oprs_t, oprs_init

    include("opr_types.jl")
    include("curl.jl")
    include("forcing.jl")
    include("nonlin.jl")
    include("integrand.jl")
    include("oprs_init.jl")

end
