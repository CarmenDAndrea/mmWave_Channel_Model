function [ SE ] = Spectral_efficiency( H,Nr,Nt, Pt_dB, M, noise_variance,h_rx,N, W, T_symbol)
% This Function calculates the spectral efficiency, assuming gaussian 
% signaling, of a MIMO system in which we transmit M symbols simultaneously 
% on the channel with channel matrix constructed as in equation (1) in:
%
% S. Buzzi, C. D'Andrea , "On Clustered Statistical MIMO Millimeter Wave Channel Simulation",
% submitted to IEEE Wireless Communications Letters
%
% For details on system model used in this function see the Section IV in
% the article listed above.
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original articles listed above.

%% INPUT PARAMETERS

% H: 3-dimensional channel matrix 

% Nr and Nt: number of receive and transmit antennas in MIMO system

%Pt_dB: transmitt power in dBW

% M: number of symbols transmitted simultaneously on channel

% noise_variance: variance of noise 

% h_rx: receive shaping filter

% N: downsampling factor of receive shaping filter

% W: System bandwidth

% T_symbol: symbol interval

%% OUTPUT PARAMETERS

% SE: spectral efficiency in bit/s/Hz


P=size(H,3); % length of channel

% select the submatrix of H for the number of antennas in the system
H=H(1:Nr,1:Nt,:); 

%% Covariance Matrix of noise

% We consider correlated noise because at each antennas we have white
% gaussian noise, but it is filtered by shaping receive filter 

% downsampling h_rx and normalizing filter energy
h_rx=h_rx(1:N:end); 
E_h=h_rx'*h_rx;
h_rx=h_rx/sqrt(E_h);

T=numel(h_rx);
rw1=zeros(1,Nr*P);
h_rx_shift=[h_rx;zeros(T,1)];
for ll=1:T
    h_rx_shiftl=circshift(h_rx_shift,[ll-1 0]);
    h_rx_shiftl=h_rx_shiftl(1:T);
    rw1(1,ll)=noise_variance*h_rx'*conj(h_rx_shiftl);
end   
rw2=[rw1(1,1) conj(rw1(2:end))];
C_w=toeplitz(rw2,rw1);

%% Precoding and Combining

% For the choice of precoding and combining matrix we identify the matrix 
% with dimensions [Nr x Nt] with maximum Frobenius norm, between P matrices
% that compose H. So we choice as precoding matrix with dimensions [Nt x M] 
% a matrix contains the rigth eigenvectors of the channel matrix 
% with the maximum Frobenius norm and as combining matrix with 
% dimensions [Nt x M] the left eigenvectors of the channel matrix with
% maximum Frobenius norm.

H_frobenius=zeros(P,1);
for p=1:P
   H_frobenius(p,1)=norm(H(:,:,p),'fro'); 
end
[~,index_max]=max(H_frobenius);
[U,~,V]=svd(H(:,:,index_max)); 
Q_precoding=V(:,1:M);
D_postcoding=U(:,1:M);

%% Power in naturals
Pt=10.^(Pt_dB/10);

%% Construction of matrix A
A_shift=[];
for t=1:P
   A_shift=[A_shift D_postcoding'*H(:,:,t)*Q_precoding];
end
A_shift=[A_shift zeros(M,M*(P-1))];
A=zeros(M*P,M*(2*P-1));
for row=0:P-1
    A(row*M+1:(row+1)*M,:)=circshift(A_shift,[0 row*M]);
end
%% Construction of matrix B
B=kron(eye(P),D_postcoding');

%% Construction of LMMSE matrix
R_r=Pt/M*(A*A')+B*C_w*B';
I_sr=transpose([zeros(M,M*P) eye(M) zeros(M,M*(2*P-1-P-1))]);
R_rs=Pt/M*A*I_sr;
E=R_r\R_rs;
A_sign=A(:,M*P+1:M*P+M);
A_int=A;
A_int(:,M*P+1:M*P+M)=[];
R=E'*(B*C_w*B'+Pt/M*(A_int*A_int'))*E;
%% Calculation of Rate with Gaussian symbols (in bit per channel use)
Rate= log2(det(eye(M)+pinv(R)*(Pt/M*E'*(A_sign*A_sign')*E)));

%% Calculation of Spectral Efficiency
SE=real(Rate/(W*T_symbol));
end

