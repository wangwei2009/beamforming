%   An example MatLAB code to compute deconvolution of signals recorded at
%   the surface and six level downhole array. Downhole signals were
%   deconvolved with the signal at the surface in frequency domain.
%
%   Devonvolution function is available at 
%   https://www.mathworks.com/matlabcentral/fileexchange/60644-deconvolution-of-two-discrete-time-signals-in-frequency-domain
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS 
%   OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
%   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
%   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE   
%   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
%   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
%   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
%   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
%   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
%   $Revision: 1.0 $  $Date: 2016/09/01 08:00:00 $
clc; close all; clear all
load test.mat; % Input data containing signals from seven accelerometers 
% Low pass filter input data
[B,A] = butter(5,10/100,'low');
HNN = filter(B,A,HNN);
% Depth of borehole levels in meter
depth = [61, 45.4, 30.5, 18.3, 10.7, 4.6];
% Plotting
figure(1);
fsz = 17; % Font size
set(gcf,'position',[300 283 400 400]);
set(gca,'TickLength',[.01 .01]);
hold on;
nlayers = 5; % Data at -4.6 m  depth will be ignored, which reduces to number of layers to five from six
% Resample data to higher sampling rate for improved pick peaking 
HNN = resample(HNN,5,1);
dt_org = 0.005; % delta t of 0.005 s (200 samples-per-second)
dt = dt_org/5;  % new delta t, which is 1000 samples-per-second
u2 = HNN(20001:80000,1);  % use windowed surface data for deconlution 
for j = 1:nlayers
    u1 = HNN(20001:80000,8-j);
    [t,w] = deconvolution(u1,u2);   
    wm = max(abs(w));
    w = w/wm; % normalize deconvolved wave amplitude
    t = t*dt;
    w = w - depth(j)/10;    
    plot(t,w,'k','LineWidth',1); % grid on;
end
xlim([-0.3 0.3]); ylim([-6.5 0]); box on;
set(gca,'Xtick',[-0.3:0.1:0.3]);
set(gca,'XTickLabel',{'-0.3','-0.2','-0.1','0','0.1','0.2','0.3'},'FontSize',[fsz],'fontname','times');
set(gca,'YTickLabel',{'-60','','-40','','-20','','Surface'},'FontSize',[fsz],'fontname','times');
xlabel('Time, s','FontSize',fsz,'FontWeight','normal','fontname','times');
ylabel('Depth, m','FontSize',fsz,'FontWeight','normal');
set(gca,'FontSize',fsz,'fontname','times');