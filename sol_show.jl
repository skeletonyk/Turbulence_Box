module sol_show
    using Plots
    export sshow, sol_one

    gr(show = true)
    function sshow(sols, n)

        w = sols.sol[n]
        N = size(w,3)
        ω = real(irfft(w, N, [1,2,3]))
        w_norm = sqrt.(real.(ω[:,:,:,1]).^2 + real.(ω[:,:,:,2]).^2+ real.(ω[:,:,:,3]).^2)
        z = w_norm[:,:,Int(size(w,3)/2)]
        display(contour(z, fillrange=true))
    end
    function sol_one(w)
        N = size(w,3)
        ω = real(irfft(w, N, [1,2,3]))
        w_z = real.(ω[:,:,size(w,3)>>1,3])
        display(contour(w_z, fillrange=true))

    end
end
