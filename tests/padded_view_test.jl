using PaddedViews
using BenchmarkTools
function pad_test(a, N, b, c, f)
    N_4 = Int(N /4 )
    println("----------------")
    @time b = PaddedView(0.0, a, (N + N_4*2, N + N_4*2, N + N_4*2), (N_4+1, N_4+1, N_4+1))
    d .= b
    @time c .=   f_plan * b
    #A_mul_B!(c, f, b)
    nothing
end

N = 128
N_4 = Int(N /4)
a = rand(N,N,N)
b = PaddedView(0.0, a, (N + N_4*2, N + N_4*2, N + N_4*2), (N_4+1, N_4+1, N_4+1))
f_plan = plan_fft(b, flags = FFTW.MEASURE)
#A_mul_B!(c, f_plan, b)
ran = rand(N_4*6, N_4*6, N_4*6)
@btime f_plan * b
@btime f_plan *ran

#@btime pad_test(a,N,b,c, f_plan)
