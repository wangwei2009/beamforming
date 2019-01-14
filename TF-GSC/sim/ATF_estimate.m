function [ h_est] = ATF_estimate(z)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明

M = size(z,2);
N_FFT = 256;
h_est = zeros(N_FFT,M-1);
segmentNum = 13;
segmentLen = 1000;
frameLen = 1024;
N = frameLen/4;

L = length(h_est(:,1));
for m = 2:M
    % matlab built-in transfer function estimator H = Pxy/Pxx
    txy =tfestimate(z(:,1),z(:,m),N_FFT,[],[],'twoside');
    
    h_est(:,m-1) = ifft(txy);
        
    % non-casual
    h_est(:,m-1) = fftshift(h_est(:,m-1));
    shift = length(h_est(:,m-1))/2;
    h_est(:,m-1) = h_est((L/2+1)-shift:(L/2+1)+shift-1,m-1);
end

