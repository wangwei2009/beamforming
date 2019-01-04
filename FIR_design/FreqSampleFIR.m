%   ������λFIR,h(n)Ϊʵ���У���h(n) = ��h(N-1-n),��h(n)�Գƣ���Գƻ�ż�Գƣ�
%   NΪ������h(n)ż�Գƣ�����A��w)����ż�Գ�
%                           A(k) = A(N-k) ,k = 0,1...N-1
%                       ��λphy = -w*(N-1)/2
%                                                 = -2*pi*k/N*(N-1)/2
%                                                 = -k*pi*(N-1)/N
%   NΪż����h(n)ż�Գƣ�����A(w)������Գ�
%                           A(k) = -A(N-k) ,k = 0,1...N-1
%                       ��λphy = -k*pi*(N-1)/N
%
%   H(w) = A(w)*exp(j*phy)
%        = A(w)*exp(-j*k*pi*(N-1)/N)  ����H(w)���㹲��Գ���
%
%
%
%%
close all
N = 33;    % NΪ����
fs = 16000;
fc = 4000.0;               %cut-off frequency
wc = fix(fc*N/fs)+1;

k = 0:N-1;
Hd = [ones(1,9),0.4,zeros(1,14),0.4,ones(1,8)]; %���Ⱥ���ż�Գ�
A = exp(-1j*pi*k*(N-1)/N);
Hd1 = Hd.*A;
hd1 = ifft(Hd1);
figure,freqz(hd1)

%%
N = 32;  % NΪż��
k = 0:N-1;
wc = 8;
Hd = [ones(1,9),0.4,zeros(1,13),-0.4,-1*ones(1,8)];% ���Ⱥ�����Գ�
A = exp(-1j*pi*k*(N-1)/N);
Hd2 = Hd.*A;
hd2 = ifft(Hd2);
figure,freqz(hd2)

