% This Matlab script can be used to generate Figure 3 in:
%
% S. Buzzi, C. D'Andrea , "A Clustered Statistical MIMO Millimeter Wave
% Channel Model"
%
% In this Figure we show the CDFs of Spectral Efficiency in mmWave
% Channel with channel response generate as in expression (1) varying
% number of symbols transmit simultaneously on the MIMO channel with Nt=30
% and Nr=20 with 0 dBW as value of transmit power, with two different link lengths.
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
Zt=6; % number of transmit antennas on the z-axis of planar array
Yr=5; % number of receiver antennas on the y-axis of planar array
Zr=4; % number of receiver antennas on the z-axis of planar array

Nr=Yr*Zr; % number of antennas of the transmit array
Nt=Yt*Zt; % number of antennas of the receive array

f=28e09; % carrier frequency

% Positions of transmitter and receiver in 3-D plane
TX_pos=[0 0 7];
RX_pos_1=[10 0 1];
RX_pos_2=[60 0 1];

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

% directions of maximum in the sectors
dir_max_theta_t=0;
dir_max_theta_r=0;
dir_max_phy_t=0;
dir_max_phy_r=180;
dir_max=[dir_max_theta_t  dir_max_theta_r;dir_max_phy_t dir_max_phy_r];


%% Parameters of the filter 
% We consider RRC pulses as transmitt and receive shaping pulses 

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

%% Calculation of Spectral Efficiency
N_channels=10000;
Pt_dB=0;
% number of symbols transmitting on the channel simultaneously
M_1=2; 
M_2=4;
M_3=6;
M_4=8;

% Initializing variables
SE1_1=zeros(N_channels,1);
SE2_1=zeros(N_channels,1);
SE3_1=zeros(N_channels,1);
SE4_1=zeros(N_channels,1);
SE1_2=zeros(N_channels,1);
SE2_2=zeros(N_channels,1);
SE3_2=zeros(N_channels,1);
SE4_2=zeros(N_channels,1);


for ch=1:N_channels
    % Generation of Channel matrix LTI with RX_pos_1
    H_frequency_selective_LTI=Generate_Channel_frequency_selective_LTI(f,TX_pos,RX_pos_1,Yt,Zt,Yr,Zr,Limitsdeg_t,Limitsdeg_r,dir_max,h_r_t,Ts,Tc);

    % Spectral Efficiencies
    SE1_1(ch,1)= Spectral_efficiency( H_frequency_selective_LTI,Nr,Nt, Pt_dB, M_1 ,noise_variance,rrc_r,N,W,T_symbol);
    SE2_1(ch,1)= Spectral_efficiency( H_frequency_selective_LTI,Nr,Nt, Pt_dB, M_2 ,noise_variance,rrc_r,N,W,T_symbol);
    SE3_1(ch,1)= Spectral_efficiency( H_frequency_selective_LTI,Nr,Nt, Pt_dB, M_3 ,noise_variance,rrc_r,N,W,T_symbol);
    SE4_1(ch,1)= Spectral_efficiency( H_frequency_selective_LTI,Nr,Nt, Pt_dB, M_4 ,noise_variance,rrc_r,N,W,T_symbol);
    
     % Generation of Channel matrix LTI with RX_pos_2
    H_frequency_selective_LTI_2=Generate_Channel_frequency_selective_LTI(f,TX_pos,RX_pos_2,Yt,Zt,Yr,Zr,Limitsdeg_t,Limitsdeg_r,dir_max,h_r_t,Ts,Tc);
    
    SE1_2(ch,1)= Spectral_efficiency( H_frequency_selective_LTI_2,Nr,Nt, Pt_dB, M_1 ,noise_variance,rrc_r,N,W,T_symbol);
    SE2_2(ch,1)= Spectral_efficiency( H_frequency_selective_LTI_2,Nr,Nt, Pt_dB, M_2 ,noise_variance,rrc_r,N,W,T_symbol);
    SE3_2(ch,1)= Spectral_efficiency( H_frequency_selective_LTI_2,Nr,Nt, Pt_dB, M_3 ,noise_variance,rrc_r,N,W,T_symbol);
    SE4_2(ch,1)= Spectral_efficiency( H_frequency_selective_LTI_2,Nr,Nt, Pt_dB, M_4 ,noise_variance,rrc_r,N,W,T_symbol);
end

%% Calculation of Empirical CDFs of Spectral Efficiencies

[ x1_1,CDF1_1] = Empirical_CDF(SE1_1);
[ x2_1,CDF2_1] = Empirical_CDF(SE2_1);
[ x3_1,CDF3_1] = Empirical_CDF(SE3_1);
[ x4_1,CDF4_1] = Empirical_CDF(SE4_1);
[ x1_2,CDF1_2] = Empirical_CDF(SE1_2);
[ x2_2,CDF2_2] = Empirical_CDF(SE2_2);
[ x3_2,CDF3_2] = Empirical_CDF(SE3_2);
[ x4_2,CDF4_2] = Empirical_CDF(SE4_2);

%% Save Results

save('Results_Figure_3');

%% Figure 3
figure
plot(x1_1,CDF1_1,'LineWidth',2);
hold on
plot(x2_1,CDF2_1,'r','LineWidth',2);
hold on
plot(x3_1,CDF3_1,'g','LineWidth',2);
hold on
plot(x4_1,CDF4_1,'k','LineWidth',2);
hold on
plot(x1_2,CDF1_2,'b--','LineWidth',2);
hold on
plot(x2_2,CDF2_2,'r--','LineWidth',2);
hold on
plot(x3_2,CDF3_2,'g--','LineWidth',2);
hold on
plot(x4_2,CDF4_2,'k--','LineWidth',2);
grid on
xlabel('Spectral efficiency [bit/s/Hz]','FontSize',12);
ylabel('CDF','FontSize',12);
AX=legend('M=2,d=10 m','M=4,d=10 m','M=6,d=10 m','M=8,d=10 m','M=2,d=60 m','M=4,d=60 m','M=6,d=60 m','M=8,d=60 m','Location', 'SouthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',12)
lim_x=max([max(x1_1) max(x2_1) max(x3_1) max(x4_1) max(x1_2) max(x2_2) max(x3_2) max(x4_2)]);
xlim([0 lim_x]);
toc
