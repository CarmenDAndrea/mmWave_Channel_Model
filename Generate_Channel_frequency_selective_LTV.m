function H=Generate_Channel_frequency_selective_LTV(f,TX_pos,RX_pos,Yt,Zt,Yr,Zr,Limitsdeg_t,Limitsdeg_r,chi_r_deg,psi_r_deg,a,Ts,Tc,v,Ncl_min,Ncl_max,Nray)
% Function for the generation of a MIMO channel matrix at mmWaves
% frequencies in case of linear time variant channel.This function
% refers to equation (11) in:
%
% S. Buzzi, C. D'Andrea , "A Clustered Statistical MIMO Millimeter Wave
% Channel Model", submitted to IEEE Wireless Communications Letters
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original article listed above.

%% INPUT PARAMETERS
% f: Carrier Frequency (e.g., 28e9 Hz)

% TX_pos: 3-dimensional vector for the position of the transmitter (e.g.,
% [0 0 7] for a transmitter placed at the origin at 7 meters height)

% RX_pos: 3-dimensional vector for the position of the receiver (e.g.,
% [30 0 1] for a transmitter at distance 30m from the transmitter at 
% 1 meters height)

% Yt and Zt: number of antennas on the y and x axes of the transmitter
% planar array (e.g., 10 and 5 for a planar array with 50 antennas)

% Yr and Zr: number of antennas on the y and x axes of the receiver
% planar array (e.g., 10 and 8 for a planar array with 80 antennas)

% Nray: constant number of rays for each cluster( e.g., Nray=8)

% Ncl_min: minimum number of clusters(e.g, Ncl_min=10)

% Ncl_max: maximum number of clusters(e.g. Ncl_max=50)

% Limitsdeg_t: 2X2 matrix contains maximum and minimum angles ,in degrees,
% in elevation and azimuth for the transmitter (e.g., [45 -45;90 -90])

% Limitsdeg_t: 2X2 matrix contains maximum and minimum angles ,in degrees,
% in elevation and azimuth for the receiver (e.g., [45 -45;90 -90])

% chi_r: orientation of receiver in elevation

% psi_r: orientation of receiver in elevation

% a: is a vector contains the samples of the convolution between the
% transmit and receive shaping pulses sampling with step Ts (e.g., raised
% cosine pulse, for the convolution between two RRC shaping pulse)

% Ts: sampling time for the filter 

% Tc: sampling time for the channel( e.g, 1e-8 s)

% v: is the constant speed of receiver in km/h (e.g., 100 km/h) 

% N.B: Tc should never be smaller than Ts to avoid bad precision results. An
% interpolation routine is indeed used in order to come with missing
% samples


%% OUTPUT Parameter
% H: 4-dimensional array of dimension N_R x N_T x P x length(tt) representing the
% samples of the matrix-valued channel impulse response taken with a 
% sampling frequency equal to 1/Tc (e.g., Tc=1e-8 s) varying the time tt in
% the interval [0,max(tau)] with step Tc;

%% Initial configurations

Nr=Yr*Zr; % number of antennas of the receiver array
Nt=Yt*Zt; % number of antennas of the transmitter array

lambda=3e08/f; % wavelength
d_array=lambda/2; % distance between elements of planar array( we consider half-wavelength arrays)
k=2*pi/lambda; % wavenumber
kd=k*d_array;

path_loss_exponent=3.3; % path loss exponent in NLOS condition ( alpha in article listed above)
beta=50;% path los at 1 m in dB

% standard deviations of Laplacian distribution ( we consider 5°) 
dev_standard_phit=5/180*pi;
dev_standard_thetat=5/180*pi;

% variance of the complex gain of the paths
alpha_variance=1;

% heigths of receiver and transmitter
hr=RX_pos(3);
ht=TX_pos(3);

% coordinate in the 3D space with TX in the origin of reference system

x=RX_pos(1);
y=RX_pos(2);
z=RX_pos(3)-TX_pos(3);

% Polar Coordinate

d=sqrt(x^2+y^2+z^2); % distance between transmitter and receiver
theta_r=atan(z/d);
phi_r=atan(y/x);
% distance between receiver and transmitter on the (x,y) plane
distance1=abs(TX_pos(1)-RX_pos(1));


%% Calculation of Number of scattering cluster depending from distance between transmitter and receiver

% We use the equation (2) in the article listed above

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

chi_r=chi_r_deg/180*pi;

psi_r=psi_r_deg/180*pi;
%% Inizialization of local variables and computation of environment geometrical features

tau_matrix=zeros(Ncl,Nray);
phit=zeros(Ncl,Nray);
phir=zeros(Ncl,Nray);
thetat=zeros(Ncl,Nray);
thetar=zeros(Ncl,Nray);
r_path=zeros(Ncl,Nray);

for ii=1:Ncl 
    % Generation of position in azimuth and elevation of the i-th cluster
    % as random variables with uniform distribution in appropriate
    % intervals
    phit_pos_i=rand*2*pi;  % uniform in [0,2*pi]
    thetat_pos_i=rand*pi-pi/2; % uniform in[-pi,pi]
    
    % Definition of Rmax as the maximum r-position of the cluster to avoid
    % considering a cluster underground 
    
    Rmax=sqrt(ht^2+(ht*tan(pi/2+thetat_pos_i))^2);
    
    % Definition of r-position of the cluster
    if thetat_pos_i>0
        r=rand*(7*d/4 -1)+1; % uniform random variable in (1,7/4*d)
    else        r=rand*(min(7*d/4,Rmax)-1)+1; % uniform random variable in (1,(min(7*d/4,Rmax)) to avoid considering a cluster underground
    end
    
    
    for ll=1:Nray
        
        % Generation of the position of l-th path in the i-th cluster as a
        % Laplace random variable with mean as the position of the i-th
        % cluster
        phit(ii,ll)=Laplace_distribution(dev_standard_phit)+phit_pos_i; 
        thetat(ii,ll)=Laplace_distribution(dev_standard_thetat)+thetat_pos_i;

                
        % Calculaton of the azimuth and elevation of the position of i-th
        % cluster respect the receiver using equations (4) and (5) in article
        % listed above
        thetar(ii,ll)=-atan((hr-ht-r*sin(thetat(ii,ll)))/(distance1-r*cos(thetat(ii,ll))))-psi_r;
        phir(ii,ll)=pi-atan(r*sin(phit(ii,ll))/(distance1-r*cos(phit(ii,ll))))-chi_r;
        
        % Calculation of the length of the path ii-ll as in equations (6)
        % in article listed above
        
        c1l=abs(ht-hr)+r*sin(thetat(ii,ll));
        c2l=d-r*cos(thetat(ii,ll));
        r_rxl=sqrt(c1l^2+c2l^2);
        r_path(ii,ll)=r+r_rxl;
        % calculation of delay of the path ii-ll 
        
        tau_matrix(ii,ll)=(r_path(ii,ll))/3e08;

    end
end


tt=0:Tc:max(max(tau_matrix));
dim_t=length(tt);

% calculation of the temporary complex path gain
alpha_temp=(randn(Ncl,Nray,dim_t)+1j*randn(Ncl,Nray,dim_t))*sqrt(alpha_variance/2);

% Introduction of correlation in time with a correlation matrix

R=zeros(dim_t);
for hh=1:dim_t
    for pp=1:dim_t
        R(hh,pp)=0.995^(abs(hh-pp));
    end
end
% Generation of complex path gains incorrelated for different paths but
% correlated in time for each path
alpha=zeros(Ncl,Nray,dim_t);
for ii=1:Ncl
    for ll=1:Nray
        temp=zeros(dim_t,1);
        for ff=1:dim_t
            temp(ff,1)=alpha_temp(ii,ll,ff);
        end
        alpha(ii,ll,:)=sqrtm(R)*temp;
    end 
end
%% Computation of channel matrix

PP=ceil((max(max(tau_matrix))/(Tc)) + length(a)); % definition of temporary length for the matrix H_temp
H_temp=zeros(Nr,Nt,PP,dim_t);

for t=1:dim_t
    %% Consideration of LOS path as in equation (10) in the article listed above
if d<5
    % if distance is less than 5 meters there is a LOS path with
    % probability 0.5
    u=rand;
    if u<=0.5
        psi=rand*2*pi; % uniform random variable in [0,2*pi]
        alpha_LOS=exp(1j*psi);
        attenuation_LOS=sqrt(10^(-0.1*beta)*d^(-path_loss_exponent));
        theta_t_LOS=theta_r;
        theta_r_LOS=-theta_r-psi_r;
        phi_t_LOS=phi_r;
        phi_r_LOS=-phi_r-chi_r;
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
        theta_t_LOS=theta_r;
        theta_r_LOS=-theta_r-psi_r;
        phi_t_LOS=phi_r;
        phi_r_LOS=-phi_r-chi_r;
        tau_LOS=d/3e8;
        at=Array_response(Yt,Zt,phi_t_LOS,theta_t_LOS,kd);
        ar=Array_response(Yr,Zr,phi_r_LOS,theta_r_LOS,kd);
        v_ms=v*1e3/3600; % express the speed of receiver in m/s
        ni_LOS=-f*v_ms/3e8*cos(theta_r_LOS)*cos(phi_r_LOS); 
        Q_LOS=gamma*alpha_LOS*attenuation_LOS*exp(-2j*pi*ni_LOS*tt(t))*ar*at';
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
            
            % Lambda Functions as in equation (8)  in article
            % listed above
            Lambdat_il=Function_Lambda(phit(ii,ll),thetat(ii,ll),phit_min,phit_max,thetat_min,thetat_max);
            Lambdar_il=Function_Lambda(phir(ii,ll),thetar(ii,ll),phir_min,phir_max,thetar_min,thetar_max);
            
            % Calculation of Doppler shift as in equations (12) in article
            % listed above
            v_ms=v*1e3/3600; % express the speed of receiver in m/s
            ni_il=-f*v_ms/3e8*cos(thetar(ii,ll))*cos(phir(ii,ll)); % doppler shift
            
            % Calculation of the attenuation of each path as in equation (7)
            % in article listed above
            attenuation_il=sqrt(10^(-0.1*beta)*r_path(ii,ll)^(-path_loss_exponent));
            
            % Temporary channel matrix
            QQ=gamma*alpha(ii,ll,t)*attenuation_il*exp(-2j*pi*ni_il*tt(t))*Lambdat_il*Lambdar_il*ar*at';
            
            % product between temporary channel matrix and the interpolated
            % filter a
            gg=floor(tau_matrix(ii,ll)/Tc);
            gg_up=floor((tau_matrix(ii,ll)+(length(a)-1)*Ts)/Tc);
            % Interpolation of filter a
            Ts_step=0:Ts:(length(a)-1)*Ts;
            Tc_step=(gg+1)*Tc-tau_matrix(ii,ll):Tc:gg_up*Tc-tau_matrix(ii,ll);
            a_interp=interp1(Ts_step,a,Tc_step,'spline');
            
            % Construction of channel matrix corresponding to the path i-th
            % l-th
            for ee=1:length(a_interp)
                H_temp(:,:,ee+gg,t)= H_temp(:,:,gg+ee,t) + QQ .*a_interp(ee);
            end
        end
    end
end

% Selection of effective channel matrix with the identification of indmin
% and indmax that indicate the start and end of discrete channel response,
% so we consider the relative delay as explained in article listed above
indmin=1;
indmax=PP;

% Calculation of Frobenius norm of matrices in third dimension
Frobenius_norm=zeros(PP,1);
for pp=1:PP
    for t=1:dim_t
       Frobenius_norm(pp,1)=Frobenius_norm(pp,1)+norm(H_temp(:,:,pp,t),'fro'); 
    end
end

 % We choice as indmin the firts index of matrix H_temp in which we have a
 % non zero contribute 
for pp=1:PP
   if Frobenius_norm(pp,1)~=0
       indmin=pp;
       break
   end
end

 % We choice as indmaxn the last index of matrix H_temp in which we have a
 % non zero contribute 
for qq=0:PP-1
   if Frobenius_norm(PP-qq,1)~=0 
       indmax=PP-qq;
       break
   end
end

% Selection of effective channel matrix
H=H_temp(:,:,indmin:indmax,:);
end