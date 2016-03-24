function Lambda=Function_Lambda(phi,theta,phi_min,phi_max,theta_min,theta_max)
% This function calculates Lambda with value 1 if the l-th path in the i-th
% cluster is intercepted by the trasmitter or receiver
% and 0 otherwise. Lambda is defined in equations (10) and (11) in:
%
% S. Buzzi, C. D'Andrea , "Generation of MIMO Channels at mm-Wave
% Frequencies for 5G Systems"
%
% License: This code is licensed under the GPLv2 License.If you in any way 
% use this code for research that results in publications, please cite our
% original article listed above.

%% INPUT PARAMETERS:

% phi and theta: position in
% azimuth and elevation of the considered path;

% phi_min and phi_max: limits of the sector in azimuth (e.g. 0 and 180)

% theta_min and theta_max: limits of sector in elevation (e.g. -45 and 45)


%% OUTPUT PARAMETERS

% Lambda: scalar with value 0 or 1 as explained above.



% Initialize the value of Lambda in azimuth and elevation
Lambdaphi=0;
Lambdatheta=0;

% Consideration of two configuration of phi_min and phi_max for
% calculation of Lambdaphi( to consider periodicity of angles in radiants)
if phi_min>=0 && phi_max>=0 
   if (phi<=phi_max && phi>=phi_min) 
         Lambdaphi=1;
   end
elseif phi_min<=0 && phi_max>=0
    if phi>=0
    if phi<=phi_max || (phi>=phi_min+2*pi && phi<=2*pi) 
           Lambdaphi=1;
    end
    else
        if (phi>=phi_min) 
            Lambdaphi=1;
        end
    end
end

% calculation of Lambdatheta
if (theta<=theta_max && theta>=theta_min)
       Lambdatheta=1;
end

% Calculation of Lambda as the product of Lambdaphi and Lambdatheta,
% because it has value 1 if and only if Lambdaphi and Lambdatheta has value
% 1
Lambda=Lambdaphi*Lambdatheta;
end