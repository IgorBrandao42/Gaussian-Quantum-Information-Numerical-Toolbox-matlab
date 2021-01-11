function [I, S_tot, S] = mutual_information(V_cell)
% Calculates the mutual information for a system with covariance matrix V_cell{i}
%
% PARAMETERS:
%  V_cell - Cell where each entry is a full covariance matrix for the whole system
%
% Calculates:
% I(j)     - mutual information  for the total system of the j-th covariance matrix
% S_tot(j) - von Neumann entropy for the total system of the j-th covariance matrix
% S(i, j)  - von Neumann entropy for the i-th mode    of the j-th covariance matrix

N_CM  = length( V_cell );                 % Number of covariance matrices
N_modes = length( V_cell{1} )/2;          % Number of modes in each covariance matrix

S = zeros([N_CM, N_modes]);               % Variable to store the entropy of each mode

for j=1:N_modes                           % Loop through each mode
  V = single_mode_CM(V_cell, j);          % Get the covariance matrix for only the i-th mode
  S(:, j) = von_Neumann_Entropy(V);       % von Neumann Entropy for i-th mode of each covariance matrix
end

S_tot = von_Neumann_Entropy(V_cell);      % von Neumann Entropy for the total system of each covariance matrix

I = zeros( size(V_cell) );                % Variable to store the mutual information
for j=1:N_modes                           % Calculation of the mutual information
  I = I + S(:, j);
end
I = I - S_tot;                            % Calculation of the mutual information

end