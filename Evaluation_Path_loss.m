function [L]=Evaluation_Path_loss(path_length,f,scenario,LOS)
% Function for the evaluation of the attenuation as in equation(2) in
%
% S. Buzzi, C. D'Andrea , "On Clustered Statistical MIMO Millimeter Wave Channel Simulation",
% submitted to IEEE Wireless Communications Letters
%
% The values used here are in Table I in the paper listed above.
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original article listed above.
%
%% INPUT PARAMETERS
% path_length: the length of the path of which the function calculate the
% attenuation 

% f: is the carrier frequency

% scenario: variable that contains information about the use-case scenario,
% it assumes the values:
% - scenario==1  ==> 'Open square'
% - scenario==2  ==> 'Street Canyon'
% - scenario==3  ==> 'Indoor Office'
% - scenario==4  ==> 'Shopping mall'

% LOS: variable that is '1' if transmitter and receiver are in LOS and '0'
% otherwise

%% OUTPUT PARAMETERS

% L: attenuation of the path in the indicated scenario expressed in
% naturals

if scenario==1 % 'Open square'
    if LOS==1
        n=1.85;
        sigma=4.2;
        b=0;
        f0=1e9;
    else
        n=2.89;
        sigma=7.1;
        b=0;
        f0=1e9;        
    end
elseif scenario==2 % 'Street Canyon'
    if LOS==1
        n=1.98;
        sigma=3.1;
        b=0;
        f0=1e9;        
    else
        n=3.19;
        sigma=8.2;
        b=0;
        f0=1e9;        
    end
elseif scenario==3 % 'Indoor Office'
    if LOS==1
        n=1.73;
        sigma=3.02;
        b=0;
        f0=1e9;       
    else
        n=3.19;
        sigma=8.29;
        b=0.06;
        f0=24.2e9;        
    end 
elseif scenario==4 % 'Shopping mall'
    if LOS==1
        n=1.73;
        sigma=2.01;
        b=0;
        f0=1e9;        
    else
        n=2.59;
        sigma=7.40;
        b=0.01;
        f0=39.5e9;        
    end
end
L_dB=20*log10(4*pi*f/3e8)+10*n*(1+b*((f-f0)/f0))*log10(path_length)+normrnd(0,sigma);

L=10^(-L_dB/10);

end

