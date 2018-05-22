using FFTW
N  = 8
N_pad = 2 * N
N_2 = N>>1
ind = [collect(1:N_2); collect(N_pad - N_2+1:N_pad)]

a = rand(N,N)
fa = fft(a)

fa2 = fft(a.*a)

fa[ N_2 + 1, :] = 0
fa[ :, N_2 + 1] = 0
f_pad = complex(zeros(N_pad, N_pad))

f_view = view(f_pad,ind,ind)
f_view .= fa
b1 = fft((ifft(f_pad)).*(ifft(f_pad)))
b1_cut = b1[ind, ind]
