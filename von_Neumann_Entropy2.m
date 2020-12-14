function result = von_Neumann_Entropy2(V_cell, j, k)
% Calculation of the von Neumann entropy for a bipartite system
%
% PARAMETERS:
%  t - times at which the von Neumann entropy will be evaluated
%
% MATHEMATICAL DESCRIPTION
%  Calculation of the von Neumann entropy for a bipartite system
%
% The covariance matrix for a bipartite subsystem is of the form:
%
%      |  A      C  |
%  V = |            | ,
%      |  C^T    B  |
%
% where A, B, C are 2 by 2 matrices and C^T is the transpose of matrix C.
% See covariance_matrix2.m for more deatils on the calculation.
%
% For a bipartite system, the von Neumann entropy (S{N}) is a
% function of the symplectic eigenvalues of the covariance matrix
% (\nu_{minus} and \nu_{plus}):
%
% S{N} = f(\nu_plus) + f(\nu_minus)
%
% where f(x) = [(x+1)/2]*log((x+1)/2) - [(x-1)/2]*log((x-1)/2) ,
% \nu_plus  = sqrt(\delta/2 + sqrt(\delta^2 - 4det(V)/2))    and
% \nu_minus = sqrt(\delta/2 + sqrt(\delta^2 - 4det(V)/2))      ,
%
% where \Delta = det(A) + det(B) + 2det(C)


% Allocate space in memory to save the value of the entropy at each time
Neuman_entropy = zeros([length(V_cell), 1]);

% Change from my label to the index
j = 2*j+1;
k = 2*k+1;

% Parallel loop to save time
parfor i=1:length(V_cell)
  
  V = V_cell{i}; % Take the full Covariance matrix
  
  A = V(j:j+1, j:j+1); % Only look for the submatrix corresponding to the desired subsystem
  B = V(k:k+1, k:k+1);
  C = V(j:j+1, k:k+1);
  V = [A, C; C.', B];
  
  % Auxiliar variable
  Delta = det(A) + det(B) + 2.0*det(C);
  
  % Smallest of the symplectic eigenvalues of the partially transposed covariance matrix
  n_minus = sqrt( Delta/2.0 - sqrt( Delta^2 - 4.0*det(V) )/2.0 );
  
  % Biggest of the symplectic eigenvalues of the partially transposed covariance matrix
  n_plus  = sqrt( Delta/2.0 + sqrt( Delta^2 - 4.0*det(V) )/2.0 );
  
  % Calculate the entropy at the current time
  Neuman_entropy(i) = func(n_plus) + func(n_minus);
  
end

% Return the entropies
result = Neuman_entropy;

end


