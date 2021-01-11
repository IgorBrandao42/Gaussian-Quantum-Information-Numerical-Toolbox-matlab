function lambda = symplectic_eigenvalues(V, Omega)
% Calculates the sympletic eigenvalues of a covariance matrix V with symplectic form Omega
% Finds the absolute values ofthe eigenvalues of i\Omega V and removes repeated entries
%  
% PARAMETERS:
%   V - covariance matrx whose symplatic eigenvalues are to be found
%   Omega (optional) - symplectic form, if this parameter is not passed, it finds the standard symplectic form

N = size(V, 1)/2;                        % Number of symplectic eigenvalues (modes)

if nargin == 1
  omega = [[0, 1]; [-1, 0]];             % Auxiliar variable
  Omega = [];                            % Symplectic form matrix
  for i=1:N                              % Build the symplectic form
    Omega = blkdiag(Omega, omega);
  end
end

H = 1i*Omega*V;                          % Auxilir matrix
lambda_0 = abs( eig(H) );                % Absolute value of the eigenvalues of the previous matrix

lambda = zeros(N, 1);                    % Variable to store the symplectic eigenvalues

for i=1:N                                % Loop over the non-repeated entries of lambda_0
  lambda(i) = lambda_0(1);               % Get the first value on the repeated array
  lambda_0(1) = [];                      % Delete it
  
  [~, idx] = min( abs(lambda_0-lambda(i)) ); % Find the next closest value on the array (repeated entry)
  lambda_0(idx) = [];                    % Delete it too
end

end

% N=2 analytical expressions for the symplectic eigenvalues
% A = V(1:2, 1:2);                       % Make use of its submatrices
% B = V(3:4, 3:4);
% C = V(1:2, 3:4);
%   
% Delta = det(A) + det(B) + 2.0*det(C);  % Auxiliar variable
%
% n_minus = sqrt( Delta/2.0 - sqrt( Delta^2 - 4.0*det(V) )/2.0 ); % Smallest of the symplectic eigenvalues of the partially transposed covariance matrix
% 
% n_plus  = sqrt( Delta/2.0 + sqrt( Delta^2 - 4.0*det(V) )/2.0 ); % Biggest of the symplectic eigenvalues of the partially transposed covariance matrix