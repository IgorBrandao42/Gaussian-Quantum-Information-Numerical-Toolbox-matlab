function W = wigner(x_mean, V, x)
% Calculates the wigner function for a gaussian quantum state with mean x_mean and
% covariance matrix V at phase space points x. See Rev. Mod. Phys. 84, 621.
%
% PARAMETERS:
%   x_mean - array with mean value of the quadratures in the gaussian state
%   V      - covariance matrix of the gaussian state
%   x      - array with points in phase space, x(:, i) is the i-th point

if size(x_mean, 2) == length(V)          % Just make sure there will be no problem with the dimensions !
  x_mean = x_mean.';
end
if size(x, 2) == length(V)
  x = x.';
end

N = length(V)/2;                         % Number of modes
W = zeros(size(x, 2), 1);                % Variable to store the calculated wigner function

for i=1:size(x, 2)
  dx = x(:, i) - x_mean;                 % x_mean(:,i) is the i-th point in phase space
  
  W_num = exp( -(dx.')*(V\dx)/2 );       % Numerator
  
  W_den = (2*pi)^N * sqrt(det(V));       % Denominator
  
  W(i) = W_num/W_den;                    % Wigner function
end

end