function [ U ] = BM_output( z,h )
%UNTITLED16 此处显示有关此函数的摘要
%   此处显示详细说明
M = size(z,2);
L = size(z,1);
Lh = size(h,1);
U = zeros(L+Lh-1,M-1);
if(mod(Lh,2)==0)
    delay = Lh/2;
else
    delay = (Lh-1)/2;
end
for i = 1:M-1
    U(:,i) = conv(z(:,1),h(:,i));
end
U = z(:,2:end) - U(1+delay:L+delay,:);


end

