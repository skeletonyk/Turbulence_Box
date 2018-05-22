w = sol.sol[1]
ω = ifft(w, [1,2,3])
k = solver_.opr.cache.con.k
con = solver_.opr.cache.con
a  = similar(w)
oprs.curl(w, a, k)
a .*= con.l_inv
b = ifft(a, [1,2,3])


con.f_view .= w
ω_pad = real(ifft(con.f_pad, [1,2,3]))

z1 = real(ω[:,:,Int(size(w,3)/2),3])
display(contour(z1, fillrange=true))


z2 = ω_pad[:,:,Int(size(ω_pad,3)/2),3]
display(contour(z2, fillrange=true))
