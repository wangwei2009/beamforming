function [ h ] = RIR_generator( s,beta)
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
d = 0.04;
r = [2   1.5 2 ;
    2-d 1.5 2;
    2-d*2 1.5 2 ;
    2-d*3 1.5 2;];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
%s = [1.94 3.5 2];              % Source position [x y z] (m)
L = [5 4 6];                % Room dimensions [x y z] (m)
%beta = 0.2;                 % Reverberation time (s)
n = beta*fs;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter

h = rir_generator(c, fs, r, s, L, beta,n, mtype, order, dim, orientation, hp_filter);

end

