% This Matlab script can be used to generate Figure 4 in:
%
% S. Buzzi, C. D'Andrea , "A Clustered Statistical MIMO Millimeter Wave
% Channel Model", submitted to IEEE Wireless Communications Letters
%
% In this Figure we show the CDFs of Spectral Efficiency in mmWave
% Channel with channel response generate as in expression (1) varying
% number of transmit and receive antennas, transmitting 4 symbols 
% simultaneously on the MIMO channel and 0 dBW transmit power.
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original article listed above.


clear all
close all
clc
tic

%%  Parameter for transmitter and receiver planar arrays 

Yt=5; % number of transmit antennas on the y-axis of planar array
Zt=8; % number of transmit antennas on the z-axis of planar array
Yr=5; % number of receiver antennas on the y-axis of planar array
Zr=8; % number of receiver antennas on the z-axis of planar array

Nr=Yr*Zr;
Nt=Yt*Zt;

f=28e09; % carrier frequency

% Positions of transmitter and receiver in 3-D plane
TX_pos=[0 0 7];
RX_pos=[30 0 1];

%% Geometrical parameters for the sectors 

% maximum and minimum angles ,in degrees, in elevation and azimuth for the
% transmitter 
thetatmax_deg=45;
thetatmin_deg=-45;
phitmax_deg=90;
phitmin_deg=-90;
Limitsdeg_t=[thetatmax_deg thetatmin_deg;phitmax_deg phitmin_deg];

% maximum and minimum angles ,in degrees, in elevation and azimuth for the
% receiver 
thetarmax_deg=45;
thetarmin_deg=-45;
phirmax_deg=90;
phirmin_deg=-90; 
Limitsdeg_r=[thetarmax_deg thetarmin_deg;phirmax_deg phirmin_deg];

% orientation of receiver in azimuth and elevation in degree
chi_r_deg=5;
psi_r_deg=5;

Nray=8; % constant number of rays for each cluster

Ncl_min=10;% minimum number of clusters

Ncl_max=50; % maximum number of cluster

%% Parameters of the filter 
% ( We consider RRC pulses as transmitt and receive shaping pulses )

R=0.22; % roll off factor
N=32; % downsampling factor
W=500e6;
T_symbol=(1+R)/W; % symbol time
Ts=T_symbol/N; % sampling time for the filter
tt=linspace(-4,4,8*N+1); 
Tc=T_symbol; % sampling time for the output of the receive filter 

%% Generation of transmit and receive shaping pulse
% We use RRC shaping filters

% RRC transmitter shaping filter
len=length(tt);
rrc_t=zeros(len,1);
for i=1:len
   t=tt(i);      
   if(t==0)
       rrc_t(i)= ( 1-R+4*R/pi ) ;
   elseif(abs(abs(t)-1/4/R)<1e-3)
       rrc_t(i)= (  cos(pi*t*(1-R))*pi*(1-R) + 4*R*cos(pi*t*(1+R)) - 4*R*t*sin(pi*t*(1+R))*(pi*(1+R))  )/(pi)/(1-3*(4*R*t)^2   ) ;   
   else
       rrc_t(i)= ( sin(pi*t*(1-R))+4*R*t*cos(pi*t*(1+R)) ) / (  pi*t*(1- (4*R*t)^2)  );
   end
end

% Normalization of filter as unitary energy filter
E_rrc_t=rrc_t'*rrc_t;
rrc_t=rrc_t/sqrt(E_rrc_t);

% RRC receive shaping filter
rrc_r=rrc_t; 

% Convolution between receiver and transmitter shaping filters and
% normalization
h_r_t=conv(rrc_t,rrc_r);
E_h=h_r_t'*h_r_t;
h_r_t=h_r_t/sqrt(E_h);

%% Noise variance

noise_figure=3; % noise figure of the receiver in dB
N0=-174;% PSD noise in dBm/Hz
noise_variance=W*10^(0.1*noise_figure)*10^(0.1*N0)*10^-3; % F*N0*B

%% Calculation of Spectral Efficiencies

N_channels=10000;
Pt_dB=0; % Transmit power in dBW

M=4;% number of symbols transmitting on the channel simultaneously

% Numbero of trasmit and receive antennas
Nr1=10;
Nt1=10;
Nr2=20;
Nt2=20;
Nr3=20;
Nt3=30;
Nr4=20;
Nt4=40;

% Initializing variables
SE1=zeros(N_channels,1);
SE2=zeros(N_channels,1);
SE3=zeros(N_channels,1);
SE4=zeros(N_channels,1);


for ch=1:N_channels
    % Generation of Channel matrix LTI
    H_frequency_selective_LTI=Generate_Channel_frequency_selective_LTI(f,TX_pos,RX_pos,Yt,Zt,Yr,Zr,Limitsdeg_t,Limitsdeg_r,chi_r_deg,psi_r_deg,h_r_t,Ts,Tc,Ncl_min,Ncl_max,Nray);
    P=size(H_frequency_selective_LTI,3);
    while isequal(H_frequency_selective_LTI,zeros(Nr,Nt,P))
        H_frequency_selective_LTI=Generate_Channel_frequency_selective_LTI(f,TX_pos,RX_pos,Yt,Zt,Yr,Zr,Limitsdeg_t,Limitsdeg_r,chi_r_deg,psi_r_deg,h_r_t,Ts,Tc,Ncl_min,Ncl_max,Nray);
        P=size(H_frequency_selective_LTI,3);
    end
    % Calculation of Spectral Efficiencies
    SE1(ch,1)= Spectral_efficiency( H_frequency_selective_LTI,Nr1,Nt1,Pt_dB, M ,noise_variance, rrc_r,N,W,T_symbol);
    SE2(ch,1)= Spectral_efficiency( H_frequency_selective_LTI,Nr2,Nt2,Pt_dB, M ,noise_variance, rrc_r,N,W,T_symbol);
    SE3(ch,1)= Spectral_efficiency( H_frequency_selective_LTI,Nr3,Nt3,Pt_dB, M ,noise_variance, rrc_r,N,W,T_symbol);
    SE4(ch,1)= Spectral_efficiency( H_frequency_selective_LTI,Nr4,Nt4,Pt_dB, M ,noise_variance, rrc_r,N,W,T_symbol);
end

%% Calculation of Empirical CDFs of Spectral Efficiencies
[ x1,CDF1] = Empirical_CDF(SE1);
[ x2,CDF2] = Empirical_CDF(SE2);
[ x3,CDF3] = Empirical_CDF(SE3);
[ x4,CDF4] = Empirical_CDF(SE4);

%% Save Results

save('Results_Figure_4');

%% Figure 3
figure
plot(x1,CDF1,'LineWidth',2);
hold on
plot(x2,CDF2,'r','LineWidth',2);
hold on
plot(x3,CDF3,'g','LineWidth',2);
hold on
plot(x4,CDF4,'k','LineWidth',2);
grid on
xlabel('Spectral Efficiency [bit/s/Hz]','FontSize',12);
ylabel('CDF','FontSize',12);
AX=legend('Nr x Nt=10 x 10','Nr x Nt=20 x 20','Nr x Nt=20 x 30','Nr x Nt=20 x 40','Location', 'SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12)
lim_x=max([max(x1) max(x2) max(x3) max(x4)]);
xlim([0 lim_x]);
toc
