function V = single_mode_CM(V_cell, k)
V = cell( length(V_cell) );              % Cell of the same size, but each entry is a single mode CM

k = 2*k+1;                               % Change from my label to the index

parfor i=1:length(V_cell)
  V_temp = V_cell{i};                    % Take the full Covariance matrix at the i-th entry
  
  V{i} = V_temp(k:k+1, k:k+1);           % Only look for the submatrix corresponding to the desired mode
end

end