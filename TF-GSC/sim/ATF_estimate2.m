function [ h_est,err ,yfbf_x12] = ATF_estimate2(z1,z2  )
%UNTITLED11 此处显示有关此函数的摘要
%   此处显示详细说明
M = 1024;
N = M/4;
z1_framed = enframe(z1, hamming(M), N,'r');
z2_framed = enframe(z2, hamming(M), N,'r');
K = size(z1_framed,1);
h_len = [];
Pz1z1 = cpsd(z1_framed',z1_framed', h_len,[],[],[],'twoside');
Pz2z1 = cpsd(z2_framed',z1_framed', h_len,[],[],[],'twoside');

b1 = sum(Pz1z1.*Pz2z1,2)/K;
b2 = sum(Pz1z1,2)/K.*sum(Pz2z1,2)/K;
a1 = sum(Pz1z1.^2,2)/K;
a2 = (sum(Pz2z1,2)/K).^2;
H = (b1-b2)./(a1-a2);
h_est = ifft(H);
L = length(h_est);
h_est = fftshift(h_est);
shift = length(h_est)/2;
h_est = h_est((L/2+1)-shift:(L/2+1)+shift-1);
% figure,plot(h_est)

%% fixed beamfomer
H = [ones(length(H),1),H];
H_mod = sum(abs(H).^2,2);
W = ifft(H./H_mod);
W = fftshift(W);
shift_w = size(W,1)/2;
yfbf_x12 = conv(W(:,1),z1)+conv(W(:,2),z2);
yfbf_x12 = yfbf_x12(1+shift_w:length(z1)+shift_w);

x2_2 = conv(z1,h_est);
%non-causal to causal,N/2 delay
x2_2 = x2_2(1+shift:length(z1)+shift);
% figure,subplot(2,1,1),plot(z2)
% hold on,plot(x2_2(1:length(z2)))
err = z2-x2_2(1:length(z2));
max(err)
% subplot(2,1,2),plot(err)
% sound(err,fs)

end

