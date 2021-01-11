function W = wigner(x_mean, V, x)
% Calculates the wigner function for a gaussian quantum state with mean x_mean and 
% covariance matrix V at phase space point x. See Rev. Mod. Phys. 84, 621.
% 
% PARAMETERS:
%   x_mean - mean value of the quadratures in the gaussian state
%   V      - covariance matrix of the gaussian state
%   x      - point in phase space

if size(x_mean, 2) > 1                   % Just make sure there will be no problem with the dimensions !
  x_mean = x_mean.';
end
if size(x, 2) > 1                        % Just make sure there will be no problem with the dimensions !
  x = x.';
end

N = length(V)/2;                         % Number of modes  

dx = x - x_mean;                         % Auxiliar variable

W_num = exp( -(dx.')*(V\dx)/2 );         % Numerator

W_den = (2*pi)^N * sqrt(det(V));         % Denominator

W = W_num/W_den;                         % Wigner function

end