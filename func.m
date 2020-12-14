function result = func(x)
% Auxiliar mathematical function
%
% MATHEMATICAL DESCRIPTION
% func(x) = (x+1)/2 * log( (x+1)/2 ) - (x-1)/2 * log( (x-1)/2 )

temp_plus  = (x + 1)/2.0;

temp_minus = (x - 1)/2.0;

result = temp_plus.*log(temp_plus) - temp_minus.*log(temp_minus);

end