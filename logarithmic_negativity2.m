function result = logarithmic_negativity2(V_cell)
% Calculation of the logarithmic negativity for a bipartite system
%
% PARAMETERS:
%    V_cell - cell where each entry is a bipartite covariance matrix
%
% MATHEMATICAL DESCRIPTION
% The covariance matrix for a bipartite subsystem is of the form:
%
%      |  A      C  |
%  V = |            | ,
%      |  C^T    B  |
%
% where A, B, C are 2 by 2 matrices and C^T is the transpose of matrix C.
%
% For a bipartite system, the logarithmic negativity (E_{N}) is a function of the
% smallest of the symplectic eigenvalues of the partially transposed covariance matrix (\tilde{\nu}_{minus}):
%
% E_{N} = max[0,-log( \tilde{\nu}_{minus} )]
%
% where \tilde{\nu}_{minus} = \sqrt( \sigma/2.0 - \sqrt( \sigma^2 - 4.0*\det(V) )/2.0 ) ,
% and   \sigma = \det(A) + \det(B) - 2\det(C)

Neg = zeros([length(V_cell),1]);         % Allocate space in memory to save the value of the logarithmic entropy for each CM

parfor i=1:length(V_cell)                % Parallel loop through every CM in the cell
  
  V = V_cell{i};                         % Take the full Covariance matrix
  
  A = V(1:2, 1:2);                       % Make use of its submatrices
  B = V(3:4, 3:4);
  C = V(1:2, 3:4);
  
  sigma = det(A) + det(B) - 2.0*det(C);  % Auxiliar variable
  
  ni = sigma/2.0 - sqrt( sigma^2 - 4.0*det(V) )/2.0 ; % Square of the smallest of the symplectic eigenvalues of the partially transposed covariance matrix
  
  if ni < 0.0                            % Manually perform a maximum to save computational time (calculation of a sqrt can take too much time and deal with residual numeric imaginary parts)
    Neg_at_t = 0.0;
  else
    ni = sqrt( real(ni) );               % Smallest of the symplectic eigenvalues of the partially transposed covariance matrix
    
    Neg_at_t = max([0, -log(ni)]);       % Calculate the logarithmic negativity at each time
  end
  
  Neg(i) = Neg_at_t;                     % Store the calculated logarithmic negativity properly
  
end

result = Neg;

end


