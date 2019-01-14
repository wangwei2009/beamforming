close all
%% generate simulation signal
source_pos = [1.94 3.5 2];              % Source position [x y z] (m)
noise_pos = [4 1.5 2];              % Noise position [x y z] (m)

h = RIR_generator( source_pos,0.2);
h = h*10;
h_n = RIR_generator( noise_pos,0.2);
h_n = h_n*1;

[s,fs] = audioread('cleanspeech.wav');
rng default  %initialize random number generator
noise = 0.5*randn(size(s));
% noise = zeros(size(s));

x1 = conv(s,h(1,:));
x2 = conv(s,h(2,:));
x3 = conv(s,h(3,:));
x4 = conv(s,h(4,:));
x = [x1,x2,x3,x4];

n1 = conv(noise,h_n(1,:));
n2 = conv(noise,h_n(2,:));
n3 = conv(noise,h_n(3,:));
n4 = conv(noise,h_n(4,:));
n = [n1,n2,n3,n4];

z = x+n;
z1 = z(:,1);
z2 = z(:,2);
z3 = z(:,3);
z4 = z(:,4);
M =size(z,2);
%audiowrite('output/z1.wav',z(:,1),fs)

%% estimate RTF
Start = 1;%14000;
End = length(x1);%58000;%length(x1);
[ h_rtf2,U2,yfbf_x12] = ATF_estimate( x1(Start:End),x2(Start:End) );
[ h_rtf3,U3,yfbf_x13] = ATF_estimate( x1(Start:End),x3(Start:End) );
[ h_rtf4,U4,yfbf_x14] = ATF_estimate( x1(Start:End),x4(Start:End) );
yfbf1 = (yfbf_x12+yfbf_x13+yfbf_x14)/(M-1);

% [ h_rtf2,U2,yfbf_x12] = ATF_estimate2( z1(Start:End),z2(Start:End) );
% [ h_rtf3,U3,yfbf_x13] = ATF_estimate2( z1(Start:End),z3(Start:End) );
% [ h_rtf4,U4,yfbf_x14] = ATF_estimate2( z1(Start:End),z4(Start:End) );
% yfbf2 = (yfbf_x12+yfbf_x13+yfbf_x14)/(M-1);

h_est = [h_rtf2,h_rtf3,h_rtf4];
%% BM output
U = BM_output( z,h_est );
%% fixed beamformer output
[ H,W,yfbf ] = fixedBf(z,h_est );
%% multi-channel ANC
sysorder = 1024;
y = zeros(size(U(:,1)));
e = zeros(size(U(:,1)));
w = zeros(sysorder,M-1);
d = yfbf;

for n = sysorder : length(x)- sysorder-1
	u = U(n:-1:n-sysorder+1,:) ;
    y(n)= sum(sum(w.* u));
    e(n) = d(n) - y(n) ;
% Start with big mu for speeding the convergence then slow down to reach the correct weights
     if n < 20
         mu=2;
     else
         mu=0.5;
     end
% Use adaptive step to reach the solution faster mu = 0.95 * 2/M*r(0)
%     mu=0.00095*2./(5*(0.001+var(u)));
% alpha = 0.01     
    a = u/(0.001+sum(sum(u.*u)) );
	w = w + mu/(0.001+sum(sum(u.*u)) ) * u * e(n) ;
    w = w + mu* e(n).*a ;%broadcasting
end 




