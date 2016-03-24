function a=Array_response(Y,Z,phi,theta,kd)
% This Function generates array response of a planar uniform array as in
% equations (12) and (13) in:
%
% S. Buzzi, C. D'Andrea , "Generation of MIMO Channels at mm-Wave
% Frequencies for 5G Systems"
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original article listed above.

%% INPUT PARAMETERS:

% Y and Z: number of antennas on the y and x axes of the
% planar array (e.g., 10 and 5 for a planar array with 50 antennas);

% phi and theta: are the position in azimuth and elevation of 
% the considered path;

% kd: product of wavenumber,2*pi/lambda, and 
% the inter-element spacing of the array.

%% OUTPUT PARAMETERS

% a: array response of planar array, column vector with Y*Z elements and
% unitary norm.

A=zeros(Y,Z); % initialize a temporary matrix 
for m=1:Y
    for n=1:Z
       A(m,n)=exp(1j*kd*((m-1)*sin(phi)*sin(theta)+(n-1)*cos(theta))); % calculate the element of temporary matrix 
    end
end
a=A(:)/sqrt(Y*Z); %dispose elements in a vector
end