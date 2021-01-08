function V = single_mode_CM(V_cell, k)
% Finds the covariance submatrix for the k-th mode from the full covariance matrix
% 
% k=0: cavity
% k>0: k-th nanoparticle
%
% PARAMETERS:
%   V_cell - Cell where each entry is a full covariance matrix  for the whole system
%   k      - index with the information about which mode is to be studied


V = cell( length(V_cell), 1 );              % Cell of the same size, but each entry is a single mode CM

k = 2*k+1;                               % Change from my label to the index

parfor i=1:length(V_cell)
  V_temp = V_cell{i};                    % Take the full Covariance matrix at the i-th entry
  
  V{i} = V_temp(k:k+1, k:k+1);           % Only look for the submatrix corresponding to the desired mode
end

end