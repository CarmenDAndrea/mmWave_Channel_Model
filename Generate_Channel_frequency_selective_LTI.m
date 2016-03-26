function H=Generate_Channel_frequency_selective_LTI(f,TX_pos,RX_pos,Yt,Zt,Yr,Zr,Limitsdeg_t,Limitsdeg_r,dir_max,a,Ts,Tc)
% Function for the generation of a MIMO channel matrix at mmWaves
% frequencies in case of linear time invariant channel. This function
% refers to equation (1) in:
%
% S. Buzzi, C. D'Andrea , "A Clustered Statistical MIMO Millimeter Wave
% Channel Model"
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original article listed above.

%% INPUT PARAMETERS
% f: Carrier Frequency (e.g., 28e9 Hz)

% TX_pos: 3-dimensional vector for the position of the transmitter (e.g.,
% [0 0 7] for a transmitter placed at the origin at 7 meters height)

% RX_pos: 3-dimensional vector for the position of the receiver (e.g.,
% [30 0 1] for a transmitter placed at distance 30m from the transmitter at 
% 1 meters height)

% Yt and Zt: number of antennas on the y and x axes of the transmitter
% planar array (e.g., 10 and 5 for a planar array with 50 antennas)

% Yr and Zr: number of antennas on the y and x axes of the receiver
% planar array (e.g., 10 and 8 for a planar array with 80 antennas)

% Ncl: Number of scattering cluster in the environment between transmitter 
% and receiver (e.g., 10 clusters)

% Nray: Number of propagation paths from each cluster (e.g., 8 rays for 
% each cluster)

% Limitsdeg_t: 2X2 matrix contains maximum and minimum angles ,in degrees,
% in elevation and azimuth for the transmitter (e.g., [45 -45;90 -90])

% Limitsdeg_t: 2X2 matrix contains maximum and minimum angles ,in degrees,
% in elevation and azimuth for the receiver (e.g., [45 -45;90 -90])

% dir_max: 2X2 matrix contains directions for the maximum of radiation
% for the arrays in the system in elevation for the transmitter and receiver
% and in azimuth for the transmitter and for the receiver (e.e., [45 45; 0
% 180])

% a: is a vector contains the samples of the convolution between the
% transmit and receive shaping pulses sampling with step Ts (e.g., raised
% cosine pulse, for the convolution between two RRC shaping pulse)

% Ts: sampling time for the filter 

% Tc: sampling time for the channel( e.g, 1e-8 s)

% N.B: Tc should never be smaller than Ts to avoid bad precision results. An
% interpolation routine is indeed used in order to come with missing
% samples


%% OUTPUT Parameter
% H: 3-dimensional array of dimension N_R x N_T x P representing the
% samples of the matrix-valued channel impulse response taken with a 
% sampling frequency equal to 1/Tc (e.g., Tc=1e-8 s) 

%% Initial configurations

Nr=Yr*Zr; % number of antennas of the receiver array
Nt=Yt*Zt; % number of antennas of the transmitter array

lambda=3e08/f; % wavelength
d_array=lambda/2; % distance between elements of planar array( we consider half-wavelength arrays)
k=2*pi/lambda; % wavenumber
kd=k*d_array;

path_loss_exponent=3.3; % path loss exponent in NLOS condition (alpha in the article listed above)
beta=50;% path los at 1 m in dB

% standard deviations of Laplacian distribution ( we consider 5°) 
dev_standard_phit=5/180*pi;
dev_standard_thetat=5/180*pi;

% variance of the complex gain of the paths
alpha_variance=1;

% heigths of receiver and transmitter
hr=RX_pos(3);
ht=TX_pos(3);

% distance between receiver and transmitter on the (x,y) plane
distance1=abs(TX_pos(1)-RX_pos(1));

% distance between transmitter and receiver on the 3-D plane
d=norm(RX_pos-TX_pos);

%% Calculation of Number of scattering cluster depending from distance between transmitter and receiver

% We use the equation (2) in the article listed above

Nray=8; % constant number of rays for each cluster
Ncl_min=10; % minimum number of clusters
Ncl_max=50; % maximum number of cluster
d_max=200; 

if d<d_max
    % if d is less than d_max the number of scattering cluster is a cubic
    % function of the distance between transmitter and receiver
    
    Ncl=Ncl_min+floor(((Ncl_max-Ncl_min)/d_max^3)*d^3);
    
else
    % if d is greater than or equal d_max the number of scattering
    % cluster is constant
    Ncl=Ncl_max;
end

% definition of constant gamma
gamma=sqrt(Nr*Nt/(Ncl*Nray));



%% Conversion of angles in radiants

phit_max=Limitsdeg_t(2,1)/180*pi;
phit_min=Limitsdeg_t(2,2)/180*pi;
thetat_max=Limitsdeg_t(1,1)/180*pi;
thetat_min=Limitsdeg_t(1,2)/180*pi;
phir_max=Limitsdeg_r(2,1)/180*pi;
phir_min=Limitsdeg_r(2,2)/180*pi;
thetar_max=Limitsdeg_r(1,1)/180*pi;
thetar_min=Limitsdeg_r(1,2)/180*pi;
dir_max_phy_t=dir_max(2,1)/180*pi;
dir_max_phy_r=dir_max(2,2)/180*pi;
dir_max_theta_r=dir_max(1,2)/180*pi;
dir_max_theta_t=dir_max(1,1)/180*pi;

%% Inizialization of local variables and generation of environment geometrical features

tau_matrix=zeros(Ncl,Nray);
phit=zeros(Ncl,Nray);
phir=zeros(Ncl,Nray);
thetat=zeros(Ncl,Nray);
thetar=zeros(Ncl,Nray);
alpha=zeros(Ncl,Nray);
r_path=zeros(Ncl,Nray);

for ii=1:Ncl 
    
    % Generation of position in azimuth and elevation of the i-th cluster
    % as random variables with uniform distribution in appropriate
    % intervals
    phit_pos_i=rand*2*pi; % uniform in [0,2*pi]
    thetat_pos_i=rand*pi-pi/2; % uniform in[-pi,pi]
    
    % Definition of Rmax as the maximum r-position of the cluster to avoid
    % considering a cluster underground
    
    Rmax=sqrt(ht^2+(ht*tan(pi/2+thetat_pos_i))^2);
    
    % Definition of r-position of the cluster
    if thetat_pos_i>0
        r=rand*(7*d/4 -1)+1; % uniform random variable in (1,7/4*d)
    else
        r=rand*(min(7*d/4,Rmax)-1)+1; % uniform random variable in (1,(min(7*d/4,Rmax)) to avoid considering a cluster underground
    end
    

   
    
    for ll=1:Nray
        
        % Generation of the position of l-th path in the i-th cluster as a
        % Laplace random variable with mean as the position of the i-th
        % cluster in azimuth and elevation for receiver and transmitter
        
        phit(ii,ll)=Laplace_distribution(dev_standard_phit)+phit_pos_i; 
        thetat(ii,ll)=Laplace_distribution(dev_standard_thetat)+thetat_pos_i;
        
        % Calculaton of the azimuth and elevation of the position of i-th
        % cluster respect the receiver as in equations (4) and (5)
        
        thetar(ii,ll)=-atan((hr-ht-r*sin(thetat(ii,ll)))/(distance1-r*cos(thetat(ii,ll))));
        phir(ii,ll)=pi-atan(r*sin(phit(ii,ll))/(distance1-r*cos(phit(ii,ll))));
        
        % Calculation of the length of the path ii-ll as in equation (6) in
        % the article listed above
        
        c1l=abs(ht-hr)+r*sin(thetat(ii,ll));
        c2l=d-r*cos(thetat(ii,ll));
        r_rxl=sqrt(c1l^2+c2l^2);
        r_path(ii,ll)=r+r_rxl;
        
        % calculation of delay of the path ii-ll 
        
        tau_matrix(ii,ll)=(r_path(ii,ll))/3e08;
        
        % calculation of the complex path gain as a complex gaussian with
        % zero mean and variance alpha_variance
        
        alpha(ii,ll)=(randn+1j*randn)*sqrt(alpha_variance/2); 
        
    end
end
%% Computation of channel matrix

PP=ceil(max(max(tau_matrix))/(Tc) + length(a)*Ts/Tc); % definition of temporary length for the matrix H_temp
H_temp=zeros(Nr,Nt,PP);

%% Consideration of LOS path as in equation (10) in the article listed above
if d<5
    % if distance is less than 5 meters there is a LOS path with
    % probability 0.5
    u=rand;
    if u<=0.5
        psi=rand*2*pi; % uniform random variable in [0,2*pi]
        alpha_LOS=exp(1j*psi);
        attenuation_LOS=sqrt(10^(-0.1*beta)*d^(-path_loss_exponent));
        theta_t_LOS=asin((hr-ht)/d);
        theta_r_LOS=asin((ht-hr)/d);
        phi_t_LOS=0;
        phi_r_LOS=pi;
        tau_LOS=d/3e8;
        at=Array_response(Yt,Zt,phi_t_LOS,theta_t_LOS,kd);
        ar=Array_response(Yr,Zr,phi_r_LOS,theta_r_LOS,kd);
        Q_LOS=gamma*alpha_LOS*attenuation_LOS*ar*at';
        % product between temporary channel matrix and the interpolated
        % filter a
        gg=floor(tau_LOS/Tc); % index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)
        gg_up=floor((tau_LOS+(length(a)-1)*Ts)/Tc);% index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)+(length(a)-1)*Ts
        % Interpolation of filter a
        Ts_step=0:Ts:(length(a)-1)*Ts;
        Tc_step=(gg+1)*Tc-tau_LOS:Tc:gg_up*Tc-tau_LOS;
        a_interp=interp1(Ts_step,a,Tc_step,'spline');
        % Construction of channel matrix corresponding to the path i-th
        % l-th
        for ee=1:length(a_interp)
            H_temp(:,:,ee+gg)= H_temp(:,:,gg+ee) + Q_LOS .*a_interp(ee);
        end
 
    end

elseif (d>=5 && d<200)
    % if distance is less than 5 meters there is a LOS path with
    % probability 0.5
    u=rand;
    if u<=0.11
        psi=rand*2*pi; % uniform random variable in [0,2*pi]
        alpha_LOS=exp(1j*psi);
        attenuation_LOS=sqrt(10^(-0.1*beta)*d^(-path_loss_exponent));
        theta_t_LOS=asin((hr-ht)/d);
        theta_r_LOS=asin((ht-hr)/d);
        phi_t_LOS=0;
        phi_r_LOS=pi;
        tau_LOS=d/3e8;
        at=Array_response(Yt,Zt,phi_t_LOS,theta_t_LOS,kd);
        ar=Array_response(Yr,Zr,phi_r_LOS,theta_r_LOS,kd);
        Q_LOS=gamma*alpha_LOS*attenuation_LOS*ar*at';
        % product between temporary channel matrix and the interpolated
        % filter a
        gg=floor(tau_LOS/Tc); % index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)
        gg_up=floor((tau_LOS+(length(a)-1)*Ts)/Tc);% index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)+(length(a)-1)*Ts
        % Interpolation of filter a
        Ts_step=0:Ts:(length(a)-1)*Ts;
        Tc_step=(gg+1)*Tc-tau_LOS:Tc:gg_up*Tc-tau_LOS;
        a_interp=interp1(Ts_step,a,Tc_step,'spline');
        % Construction of channel matrix corresponding to the path i-th
        % l-th
        for ee=1:length(a_interp)
            H_temp(:,:,ee+gg)= H_temp(:,:,gg+ee) + Q_LOS .*a_interp(ee);
        end
    end
end

for ii=1:Ncl,
    for ll=1:Nray,
        % Array responses as in equation (9) in article
        % listed above
        at=Array_response(Yt,Zt,phit(ii,ll),thetat(ii,ll),kd);
        ar=Array_response(Yr,Zr,phir(ii,ll),thetar(ii,ll),kd);
        
        % Lambda Functions as in equation (8) in article
        % listed above
        Lambdat_il=Function_Lambda(phit(ii,ll),thetat(ii,ll),dir_max_phy_t-abs(phit_min),phit_max+dir_max_phy_t,dir_max_theta_t-abs(thetat_min),thetat_max+dir_max_theta_t);
        Lambdar_il=Function_Lambda(phir(ii,ll),thetar(ii,ll),dir_max_phy_r-abs(phir_min),phir_max+dir_max_phy_r,dir_max_theta_r-abs(thetar_min),thetar_max+dir_max_theta_r);
        
        % Calculation of the attenuation of each path as in equation (7)
        % in article listed above
        attenuation_il=sqrt(10^(-0.1*beta)*r_path(ii,ll)^(-path_loss_exponent));
        
        % Temporary channel matrix
        QQ=gamma*alpha(ii,ll)*attenuation_il*Lambdat_il*Lambdar_il*ar*at';
        
        % product between temporary channel matrix and the interpolated
        % filter a
        gg=floor(tau_matrix(ii,ll)/Tc); % index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)
        gg_up=floor((tau_matrix(ii,ll)+(length(a)-1)*Ts)/Tc);% index of filter "approximately" corresponding to the delay tau_matrix(ii,ll)+(length(a)-1)*Ts
        % Interpolation of filter a
        Ts_step=0:Ts:(length(a)-1)*Ts;
        Tc_step=(gg+1)*Tc-tau_matrix(ii,ll):Tc:gg_up*Tc-tau_matrix(ii,ll);
        a_interp=interp1(Ts_step,a,Tc_step,'spline');
        
        % Construction of channel matrix corresponding to the path i-th
        % l-th
        for ee=1:length(a_interp)
            H_temp(:,:,ee+gg)= H_temp(:,:,gg+ee) + QQ .*a_interp(ee); 
        end
    end
end

% Selection of effective channel matrix with the identification of indmin
% and indmax that indicate the start and end of discrete channel response,
% so we consider the relative delay as explained in article listed above
indmin=1;
indmax=PP;

 % We choice as indmin the firts index of matrix H_temp in which we have a
 % non zero contribute 
for pp=1:PP
   if ~isequal(zeros(Nr,Nt),H_temp(:,:,pp))
       indmin=pp;
       break
   end
end

 % We choice as indmaxn the last index of matrix H_temp in which we have a
 % non zero contribute 
for qq=0:PP-1
   if ~isequal(zeros(Nr,Nt),H_temp(:,:,PP-qq))
       indmax=PP-qq;
       break
   end
end

% Selection of effective channel matrix
H=H_temp(:,:,indmin:indmax);
end