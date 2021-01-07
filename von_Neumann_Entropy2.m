function result = von_Neumann_Entropy2(V_cell)
% Calculation of the von Neumann entropy for a bipartite system
%
% PARAMETERS:
%   V_cell - cell where each entry is a bipartite covariance matrix
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

Neuman_entropy = zeros([length(V_cell), 1]); % Allocate space in memory to save the value of the entropy at each time

parfor i=1:length(V_cell)                % Parallel loop through every CM in the cell
  
  V = V_cell{i};                         % Take the full Covariance matrix
  
  A = V(1:2, 1:2);                       % Make use of its submatrices
  B = V(3:4, 3:4);
  C = V(1:2, 3:4);
  
  Delta = det(A) + det(B) + 2.0*det(C);  % Auxiliar variable
  
  n_minus = sqrt( Delta/2.0 - sqrt( Delta^2 - 4.0*det(V) )/2.0 ); % Smallest of the symplectic eigenvalues of the partially transposed covariance matrix
  
  n_plus  = sqrt( Delta/2.0 + sqrt( Delta^2 - 4.0*det(V) )/2.0 ); % Biggest of the symplectic eigenvalues of the partially transposed covariance matrix
  
  Neuman_entropy(i) = func(n_plus) + func(n_minus); % Calculate the entropy at the current time
  
end

result = Neuman_entropy;

end


