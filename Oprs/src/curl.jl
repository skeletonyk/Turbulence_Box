@inline function curl(x :: TB_data_t, y :: TB_data_t, k :: TB_data_t)
    @muladd begin
    y[1] .= x[3] .* k[2] .- x[2] .* k[3]
    y[2] .= x[1] .* k[3] .- x[3] .* k[1]
    y[3] .= x[2] .* k[1] .- x[1] .* k[2]
    end
    nothing
end

@inline function cross_3d( x :: TB_data_t, y :: TB_data_t, z :: TB_data_t)
    @muladd begin
    z[1] .= x[2] .* y[3]  .- x[3] .* y[2]
    z[2] .= x[3] .* y[1]  .- x[1] .* y[3]
    z[3] .= x[1] .* y[2]  .- x[2] .* y[1]
    end
    nothing
end
