f_pad = solver_.opr.cache.con.f_pad_half
a = rand(3,3)
fa = rfft(a)
fa_copy = copy(fa)
irf = plan_irfft(fa, 3)
c = similar(a)

println(fa)
A_mul_B!(c, irf, fa)
println(fa-fa_copy)
