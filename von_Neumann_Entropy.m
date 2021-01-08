function Entropy = von_Neumann_Entropy(V_cell)
% Calculation of the von Neumann entropy for a multipartite gaussian system
%
% PARAMETERS:
%  t - times at which the von Neumann entropy will be evaluated
%
% MATHEMATICAL DESCRIPTION
% Calculation of the von Neumann entropy for a multipartite gaussian system
% given the covariance matrix that defines the state
%
% For any gaussian system, the von Neumann entropy (S_{N}) is a
% function of the symplectic eigenvalues of the covariance matrix
% (\nu_{k}):
%
% S_{N} = \sum_{k=1}^{N} g(\nu_k)
%
% where g(x) = [(x+1)/2]*log((x+1)/2) - [(x-1)/2]*log((x-1)/2) ,
% and \nu_k are the sympletic eigenvalues of V,
% i.e., modulus of the eigenvalues of i \Omega V
%
% func.m perform the calculation of g(x)

N = size(V_cell{1}, 1)/2;                % Number of modes

omega = [[0, 1]; [-1, 0]];               % Auxiliar variable
Omega = [];                              % Symplectic form matrix
for i=1:N                                % Build the symplectic form
  Omega = blkdiag(Omega, omega);
end

Entropy = zeros([length(V_cell), 1]);    % Allocate space in memory to save the value of the entropy at each time

parfor i=1:length(V_cell)                % Parallel loop to save time
  
  V = V_cell{i};                         % Take the full Covariance matrix
  
  nu = symplectic_eigenvalues(V, Omega);
  
  Entropy(i) = sum( func(nu) );          % Calculate the entropy at the current time
  
end

end


