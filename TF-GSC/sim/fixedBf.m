function [ H,W,yfbf ] = fixedBf( z,h )
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
M = size(z,2);
L = size(z,1);
Lh = size(h,1);
if(mod(Lh,2)==0)
    delay = Lh/2;
else
    delay = (Lh-1)/2;
end
% W = zeros(L+Lh-1,M);
yfbf = zeros(L+Lh-1,1);
H = fft(fftshift(h));
H = [ones(size(H,1),1),H];
H_mod = sum(abs(H).^2,2);
W = ifft(H./H_mod);
W = fftshift(W);

for i = 1:M
    yfbf = yfbf + conv(W(:,i),z(:,i));
end
% yfbf = yfbf;

yfbf = yfbf(1+delay:length(z(:,1))+delay);
% yfbf = yfbf(1:length(z(:,1)));

end

