using FFTW
N  = 8
N_pad = 2 * N
N_2 = Int(N/2)
ind = [collect(1:N_2); collect(N_pad - N_2+1:N_pad)]

a = rand(N)
fa2 = fft(a.*a)

fa = fft(a)
#fa[ N>>1 + 1] = 0
f_pad = complex(zeros(N_pad))

f_view = view(f_pad,ind)
f_view .= fa
b1 = fft(real(ifft(f_pad)).*real(ifft(f_pad)))
b1_cut = b1[ind]*2

b2 = fft((ifft(f_pad)).*(ifft(f_pad)))
b2_cut = b2[ind]*2
