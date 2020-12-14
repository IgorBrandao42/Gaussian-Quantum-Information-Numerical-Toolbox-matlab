function result = logarithmic_negativity2(V_cell, j, k)
% Calculation of the logarithmic negativity for a bipartite system
%
% PARAMETERS:
%       t - times at which the logarithmic negativity will be evaluated
%
% MATHEMATICAL DESCRIPTION
% The covariance matrix for a bipartite subsystem is of the form:
%
%      |  A      C  |
%  V = |            | ,
%      |  C^T    B  |
%
% where A, B, C are 2 by 2 matrices and C^T is the transpose of matrix C.
% See covariance_matrix2.m for more deatils on the calculation.
%
% For a bipartite system, the logarithmic negativity (E_{N}) is a function of the
% smallest of the symplectic eigenvalues of the partially transposed covariance matrix (\tilde{\nu}_{minus}):
%
% E_{N} = max[0,-log( \tilde{\nu}_{minus} )]
%
% where \tilde{\nu}_{minus} = \sqrt( \sigma/2.0 - \sqrt( \sigma^2 - 4.0*\det(V) )/2.0 ) ,
% and   \sigma = \det(A) + \det(B) - 2\det(C)

% Allocate space in memory to save the value of the logarithmic entropy at each time
Neg = zeros([length(V_cell),1]);

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
  sigma = det(A) + det(B) - 2.0*det(C);
  
  % Square of the smallest of the symplectic eigenvalues of the partially transposed covariance matrix
  ni = sigma/2.0 - sqrt( sigma^2 - 4.0*det(V) )/2.0 ;
  
  % Manually perform a maximum to save computational time (calculation of a sqrt can take too much time and deal with residual numeric imaginary parts)
  if ni < 0.0
    Neg_at_t = 0.0;
  else
    % Smallest of the symplectic eigenvalues of the partially transposed covariance matrix
    ni = sqrt( real(ni) );                   % I am using the REAL PART, I have to make sure that this is fine instead of : ni = sqrt(ni);
    
    % Calculate the logarithmic negativity at each time
    Neg_at_t = max([0, -log(ni)]);
  end
  
  Neg(i) = Neg_at_t;
  
end

result = Neg;

end


