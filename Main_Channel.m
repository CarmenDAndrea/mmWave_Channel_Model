% This Matlab script can be used to generate channel matrices
% H_frequency_selective_LTI and H_frequency_selective_LTV respectively in
% equations (1) and (11) in:
%
% S. Buzzi, C. D'Andrea , "A Clustered Statistical MIMO Millimeter Wave
% Channel Model", submitted to IEEE Wireless Communications Letters
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
Zt=4; % number of transmit antennas on the z-axis of planar array
Yr=5; % number of receiver antennas on the y-axis of planar array
Zr=2; % number of receiver antennas on the z-axis of planar array

f=73e09; % carrier frequency

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

v=100; % speed of receiver in km/h

%% Parameters of the filter 
%  We consider RRC pulses as transmitt and receive shaping pulses 

R=0.22; % roll off factor
N=32; % downsampling factor
W=500e6;
T_symbol=(1+R)/W; % symbol time
Ts=T_symbol/N; % sampling time for the filter
tt=linspace(-4,4,8*N+1); 
Tc=2e-9; % sampling time for the output of the receive filter 

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

%% Frequency selective Channel Matrix LTI case 

H_frequency_selective_LTI=Generate_Channel_frequency_selective_LTI(f,TX_pos,RX_pos,Yt,Zt,Yr,Zr,Limitsdeg_t,Limitsdeg_r,chi_r_deg,psi_r_deg,h_r_t,Ts,Tc,Ncl_min,Ncl_max,Nray);

%% Frequency selective Channel Matrix LTV case

H_frequency_selective_LTV=Generate_Channel_frequency_selective_LTV(f,TX_pos,RX_pos,Yt,Zt,Yr,Zr,Limitsdeg_t,Limitsdeg_r,chi_r_deg,psi_r_deg,h_r_t,Ts,Tc,v,Ncl_min,Ncl_max,Nray);

toc
