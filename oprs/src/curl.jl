@inline function curl_exec(x :: TB_data_t, y :: TB_data_t, k :: TB_data_t)
    y1 = view(y,:,:,:,1)
    y2 = view(y,:,:,:,2)
    y3 = view(y,:,:,:,3)

    x1 = view(x,:,:,:,1)
    x2 = view(x,:,:,:,2)
    x3 = view(x,:,:,:,3)

    k1 = view(k,:,:,:,1)
    k2 = view(k,:,:,:,2)
    k3 = view(k,:,:,:,3)

    @muladd begin
    @. y1 = x3 * k2 - x2 * k3
    @. y2 = x1 * k3 - x3 * k1
    @. y3 = x2 * k1 - x1 * k2
    end

    nothing
end

@inline function cross_3d( x :: TB_data_real_t, y :: TB_data_real_t, z :: TB_data_real_t)
    y1 = view(y,:,:,:,1)
    y2 = view(y,:,:,:,2)
    y3 = view(y,:,:,:,3)

    x1 = view(x,:,:,:,1)
    x2 = view(x,:,:,:,2)
    x3 = view(x,:,:,:,3)

    z1 = view(z,:,:,:,1)
    z2 = view(z,:,:,:,2)
    z3 = view(z,:,:,:,3)


    @. begin
    @muladd begin
    z1 .= x2 .* y3  .- x3 .* y2
    z2 .= x3 .* y1  .- x1 .* y3
    z3 .= x1 .* y2  .- x2 .* y1
    end
    end
    nothing
end
