%close all
N = 4096;
fs = 16000;
fc = 4000.0;               %cut-off frequency
wc = fix(fc*N/fs);
if(mod(N,2))
    half_bin = (N+1)/2;
else
    half_bin = N/2;
end

%% ideal low-pass filter
Hd = zeros(1,N);
Hd(1:wc) = 1;
Hd(wc+1:N/2) = 0;

Hd(N/2+1) = 0;

Hd(N/2+2:end) = fliplr(Hd(2:N/2));

%% time-domain ideal filter ,rectwin
hd = ifft(Hd);
hd = fftshift(hd);
%figure,plot(hd)
M = 33;
hd = hd(N/2-(M-1)/2:N/2+(M-1)/2);
%% compare different window
N_FFT = N*4;   %frequency-domain interpolation
f = (1:N_FFT/2+1)*fs/N_FFT;
figure,

win = rectwin(M)';
h = hd.*win;
H = fft(h,N_FFT);
plot(f,pow2db(abs(H(1:N_FFT/2+1)))),%ylim([0,1.2]);
win = hamming(M)';
h = hd.*win;
H = fft(h,N_FFT);
hold on,plot(f,pow2db(abs(H(1:N_FFT/2+1))));
win = hann(M)';
h = hd.*win;
H = fft(h,N_FFT);
hold on,plot(f,pow2db(abs(H(1:N_FFT/2+1))));
legend('rectwin','hamming','hann')

