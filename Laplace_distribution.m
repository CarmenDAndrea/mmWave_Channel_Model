function angle=Laplace_distribution(dev_standard)
% This function generates a Laplacian random variable with zero mean and
% standard deviation dev_standard. It refers to expression (7) in:
%
% S. Buzzi, C. D'Andrea , "Generation of MIMO Channels at mm-Wave
% Frequencies for 5G Systems"
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original article listed above.

%% INPUT PARAMETER

% dev_standard: standard deviation of Laplacian random variable in radiants

%% OUTPUT PARAMETER

% angle: Laplacian random variable with zero mean and standard deviation
% dev_standard

% We generate this variable from the CDF (obtained from the integral of
% expression (10)) transforming an uniform random variable, in the interval
% [0,1], with the inverse function of calculated CDF using temporary
% variables c and d defined as follows.

c=1/((sqrt(2)*dev_standard)*(1-exp(-sqrt(2)*pi/dev_standard)));
d=sqrt(2)/dev_standard;

x=rand; % uniform random variable in [0,1]

if x<=1/2
    angle=1/d*log(d/c*x+exp(-d*pi));
else
    angle=-1/d*log(1-(x-1/2)*d/c);
end
end