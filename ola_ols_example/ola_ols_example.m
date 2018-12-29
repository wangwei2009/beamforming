h = load('RIR.txt');

[s,fs] = audioread('an103-mtms-senn4.wav');

scale = 10;
h = h*scale;
x1 = conv(s,h(1,:));
x2 = conv(s,h(2,:));
x3 = conv(s,h(3,:));
x4 = conv(s,h(4,:));

y1 = x1;
y2 = x2;
y3 = x3;
y4 = x4;
%audiowrite('output/y1.wav',y1,fs);
y = [y1,y2,y3,y4];
%audiowrite('output/y_ds.wav',sum(y,2)/4,fs);
M = size(y,2);
frameLen = 2048;
N_FFT = frameLen+length(h(1,:))-1;%4095;
L = N_FFT;
inc = frameLen;
y_ola = zeros(length(y1),M);
win = rectwin(frameLen);
H = fft(h(:,1:2048)',N_FFT);
%% Frequency domain convolution OLA
for i = N_FFT:inc:length(y1)-N_FFT
    yi = s(i-frameLen+1:i,:).*win;
    Yi = fft(yi,N_FFT).*fft(h(:,1:2048)',N_FFT);
    y_ola(i-frameLen+1:i-frameLen+N_FFT,:) = y_ola(i-frameLen+1:i-frameLen+N_FFT,:)+ifft(Yi);
end
%% Frequency domain convolution OLS
y_ols = zeros(length(y1),4);
for i = N_FFT:inc:length(y1)-N_FFT  
    yi = s(i-frameLen+1:i-frameLen+N_FFT,:);%.*win;
    Yi = fft(yi,N_FFT).*fft(h(1:4,:)',N_FFT);%.*B;
    i_Yi = ifft(Yi);
    y_ols(i-frameLen+1:i) = y_ols(i-frameLen+1:i)+i_Yi(end-frameLen+1:end);
end







