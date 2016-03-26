function [ x,y] = Empirical_CDF(Realizations)
% This Function calculate the x and y coordinates for the Empirical CDF
% from the vector of Realizations in input

%% INPUT PARAMETER

% Realizations: vector contains the realization for the calculation of
% Empirical CDF

%% OUTPUT PARAMETERS

% x: vector of abscissas

% y: vector contains the probabilities that Realizations is lower than
% abscissas


x=sort(Realizations); % sort vector in input in ascendent order
s=length(x); 
y=1/s:1/s:1; % probability realizations lower than abscissas
end

