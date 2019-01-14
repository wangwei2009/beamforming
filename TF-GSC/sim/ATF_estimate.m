function [ iRTF12,err,yfbf_x12] = ATF_estimate( x1,x2 )
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
Start = 200000;
End = 320000;
%M = size(z,2);
fs = 16000;
RTF12 =tfestimate(x1,x2,[],[],[],'twoside');
% RTF12 =tfestimate(x1,x2,[],[],[],[]);
% RTF12 = [RTF12;conj(flipud(RTF12(2:end-1)))];
iRTF12 = ifft(RTF12);
L = length(iRTF12);
shift = length(iRTF12)/2;
%% fixed beamfomer
H = [ones(length(RTF12),1),RTF12];
H_mod = sum(abs(H).^2,2);
W = ifft(H./H_mod);
yfbf_x12 = conv(fftshift(W(:,1)),x1)+conv(fftshift(W(:,2)),x2);
yfbf_x12 = yfbf_x12(1+shift:length(x1)+shift);
%% BM
% iRTF12(1)=0;
% iRTF12(length(iRTF12))=0;
%to ensure stability of the TFs ratios, we do not impose them to be causal
iRTF12 = fftshift(iRTF12);   % non-casuality
% iRTF12 = fftshift(ifft(fftshift(RTF12)));


%FIR structure constrain
shift = 4096;
iRTF12 = iRTF12((L/2+1)-shift:(L/2+1)+shift); 
% iRTF12 = iRTF12.*hamming(length(iRTF12));
delay = fix((length(iRTF12)-1)/2);
% figure,plot(real(iRTF12(1:2048)))

x2_2 = conv(x1,iRTF12);
%non-causal to causal,N/2 delay
x2_2 = x2_2(1+delay:length(x1)+delay);
figure,subplot(2,1,1),plot(x2)
hold on,plot(x2_2(1:length(x1)))
err = x2-x2_2(1:length(x1));
max(err)
subplot(2,1,2),plot(err)
sound(err,fs)

end

