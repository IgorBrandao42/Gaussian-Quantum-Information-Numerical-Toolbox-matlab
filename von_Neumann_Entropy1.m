function result = von_Neumann_Entropy1(V_cell, k)
% Calculation of the von Neumann entropy for a bipartite system
%
% PARAMETERS:
%  t - times at which the von Neumann entropy will be evaluated
%
% MATHEMATICAL DESCRIPTION
% Calculation of the von Neumann entropy for single mode of a
% N-mode covariance matrix
%
% For any gaussian system, the von Neumann entropy (S_{N}) is a
% function of the symplectic eigenvalues of the covariance matrix
% (\nu_{k}):
%
% S_{N} = \sum_{k=1}^{N} g(\nu_k)
%
% where g(x) = [(x+1)/2]*log((x+1)/2) - [(x-1)/2]*log((x-1)/2) ,
% and nu are the sympletic eigenvalues of V, 
% i.e., modulus of the eigenvalues of i \Omega V
%
% func.m perform the calculation of g(x)
%
% See section II.C.1 of Gaussian Quantum Information

omega = [[0, 1]; [-1, 0]];                   % Auxiliar variable
Omega = [];                                  % Symplectic form matrix
for i=1:1                                    % Build the symplectic form ( size(V_cell{1},1)/2 )
  Omega = blkdiag(Omega, omega);
end

Neuman_entropy = zeros([length(V_cell), 1]); % Allocate space in memory to save the value of the entropy at each time

k = 2*k+1;                                   % Change from my label to the index

parfor i=1:length(V_cell)                       % Parallel loop to save time
  
  V = V_cell{i};                             % Take the full Covariance matrix
  V = V(k:k+1, k:k+1);                       % Only look for the submatrix corresponding to the desired mode
  
  H = 1i*Omega*V;                            % Calculate the sympletic eigenvalues for this covariance submatrix
  nu = abs( eig(H) );
  nu = ( nu(1) + nu(2) )/2.0;

  Neuman_entropy(i) = sum( func(nu) );       % Calculate the entropy at the current time
  
end

result = Neuman_entropy;                     % Return the calculated entropies

end


