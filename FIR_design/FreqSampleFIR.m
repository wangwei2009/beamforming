%   线性相位FIR,h(n)为实序列，且h(n) = ±h(N-1-n),即h(n)对称（奇对称或偶对称）
%   N为奇数，h(n)偶对称，幅度A（w)满足偶对称
%                           A(k) = A(N-k) ,k = 0,1...N-1
%                       相位phy = -w*(N-1)/2
%                                                 = -2*pi*k/N*(N-1)/2
%                                                 = -k*pi*(N-1)/N
%   N为偶数，h(n)偶对称，幅度A(w)满足奇对称
%                           A(k) = -A(N-k) ,k = 0,1...N-1
%                       相位phy = -k*pi*(N-1)/N
%
%   H(w) = A(w)*exp(j*phy)
%        = A(w)*exp(-j*k*pi*(N-1)/N)  最终H(w)满足共轭对称性
%
%
%
%%
close all
N = 33;    % N为奇数
fs = 16000;
fc = 4000.0;               %cut-off frequency
wc = fix(fc*N/fs)+1;

k = 0:N-1;
Hd = [ones(1,9),0.4,zeros(1,14),0.4,ones(1,8)]; %幅度函数偶对称
A = exp(-1j*pi*k*(N-1)/N);
Hd1 = Hd.*A;
hd1 = ifft(Hd1);
figure,freqz(hd1)

%%
N = 32;  % N为偶数
k = 0:N-1;
wc = 8;
Hd = [ones(1,9),0.4,zeros(1,13),-0.4,-1*ones(1,8)];% 幅度函数奇对称
A = exp(-1j*pi*k*(N-1)/N);
Hd2 = Hd.*A;
hd2 = ifft(Hd2);
figure,freqz(hd2)

