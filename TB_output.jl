module TB_output
    using WriteVTK
    using oprs

    export write_vtk_single
    const FloatType = Float32
    const vtk_filename_noext = "structured"

    function write_vtk_single(sol, fname = "null")
       outfiles = String[]
       N = size(sol,3)
       # generate vtk_grid
       xyz = zeros(FloatType, 3, N, N, N)
       for k = 1:N, j = 1:N, i = 1:N
           xyz[1, i, j, k] = i/N
           xyz[2, i, j, k] = j/N
           xyz[3, i, j, k] = k/N
       end
       vtk = vtk_grid(fname, xyz)

       # ---------------------------------------------------------------------
       # writing vorticity
       ω = real(irfft(sol, N,[1,2,3]))

       @views vort_tuple = (
                    ω[:, :, :, 1],
                    ω[:, :, :, 2],
                    ω[:, :, :, 3])

        vtk_point_data(vtk, vort_tuple, "vorticity")

        ω_norm = sqrt.(ω[:,:,:,1].^2 + ω[:,:,:,2].^2 + ω[:,:,:,3].^2)
        vtk_point_data(vtk, ω_norm, "vort_magnitude")


       # ---------------------------------------------------------------------


        # write to outfile
        #append!(outfiles, vtk_save(vtk))

        return vtk
    end

    function write_vtk_series(sols, pvd_name)
        pvd = paraview_collection(pvd_name)
        for i = 1: length(sols.sol)
            vtk = write_vtk_single(sols.sol[i], "T_"*string(i))
            vtk_save(vtk)
            collection_add_timestep(pvd, vtk, Float64(i))
        end

        vtk_save(pvd)
        return pvd
    end

end
