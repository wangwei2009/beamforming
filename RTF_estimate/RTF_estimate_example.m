close all
%% use RIR to generate reverberation output
example_3
[s,fs] = audioread('an103-mtms-senn4.wav');
scale = 10;
h = h*scale;
x1 = conv(s,h(1,:));
x2 = conv(s,h(2,:));

%% estimate ATF
txy =tfestimate(s,x1(1:length(s)),[],[],[],[],'twoside');
itxy = ifft(txy);
itxy = itxy(1:2048);
figure,
subplot(2,1,1),plot(real(itxy)),title('estimate ATF');
subplot(2,1,2),plot(h(1,:)),title('true ATF');

x1 = conv(s,h(1,:));
x1_2 = conv(s,itxy);
% x1 = x1(384:end);
figure,subplot(2,1,1),plot(x1)
hold on,plot(x1_2(1:length(x1))),legend('x1','x1_est')
subplot(2,1,2),plot(x1-x1_2(1:length(x1))),title('residual signal');
sound(x1-x1_2(1:length(x1)),fs)

%% estimate RTF 
RTF12 =tfestimate(x1,x2,16384,[],[],[],'twoside');
% RTF12 =tfestimate(x1,x2,16384,[],[],[]);
% RTF12 = [RTF12;conj(flipud(RTF12(2:end-1)))];
iRTF12 = ifft(RTF12);
%to ensure stability of the TFs ratios, we do not impose them to be causal
iRTF12 = fftshift(iRTF12);   % non-casuality
L = length(iRTF12);
shift = 4096;
%FIR structure constrain
iRTF12 = iRTF12((L/2+1)-shift:(L/2+1)+shift); 
% figure,plot(real(iRTF12(1:2048)))

x2_2 = conv(x1,iRTF12);
%non-causal to causal,N/2 delay
x2_2 = x2_2(1+shift:length(x1)+shift);
figure,subplot(2,1,1),plot(x2)
hold on,plot(x2_2(1:length(x1)))
err = x2-x2_2(1:length(x1));
max(err)
subplot(2,1,2),plot(err)
sound(err,fs)

%% estimate RTF,proposed
M = 16384;
N = M/2;

z1 = enframe(x1, hamming(M), N,'z');
z2 = enframe(x2, hamming(M), N,'z');

K = size(z1,1);
Pz1z1 = cpsd(z1',z1', M,N,[],[],'centered');
Pz2z1 = cpsd(z2',z1', M,N,[],[],'centered');

b1 = sum(Pz1z1.*Pz2z1,2)/K;
b2 = sum(Pz1z1,2)/K.*sum(Pz2z1,2)/K;
a1 = sum(Pz1z1.^2,2)/K;
a2 = (sum(Pz2z1,2)/K).^2;
H = (b1-b2)./(a1-a2);
h_est = ifft(H);
L = length(h_est);
h_est = fftshift(h_est);
shift = 4096;
h_est = h_est((L/2+1)-shift:(L/2+1)+shift);
figure,plot(h_est)

x2_2 = conv(x1,h_est);
%non-causal to causal,N/2 delay
x2_2 = x2_2(1+shift:length(x1)+shift);
figure,subplot(2,1,1),plot(x2)
hold on,plot(x2_2(1:length(x1)))
err = x2-x2_2(1:length(x1));
max(err)
subplot(2,1,2),plot(err)
sound(err,fs)
