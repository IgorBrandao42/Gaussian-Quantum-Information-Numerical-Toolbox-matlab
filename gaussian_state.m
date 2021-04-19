classdef gaussian_state < handle                                 % class definning a nanoparticle
  properties
    R                                    % Mean quadratures
    V                                    % Covariance matrix
    
  end
  
  methods
    function obj = gaussian_state(varargin)
      % Class simulating a gaussian state with mean quadratures and covariance matrix
      % The user can explicetilly pass the first two moments of a multimode gaussian state
      % or pass a name-value pair argument to choose a single mode gaussian state
      %
      % PARAMETERS:
      %   R0 - mean quadratures for gaussian state
      %   V0 - covariance matrix for gaussian state
      %
      % NAME-VALUE PAIR ARGUMENTS:
      %   "vacuum"                        -> generates vacuum   state
      %   "thermal" , occupation number   -> generates thermal  state
      %   "coherent", complex amplitude   -> generates coherent state
      %   "squeezed", squeezing parameter -> generates squeezed state
      %
      % Calculates the initial occupation number and occupation number of the associated heat bath
      
      % TO DO (?): Make sure all entries are real entries and symmetric CM satisfiyng the uncertainty principle
      
      if (isstring(varargin{1}) || ischar(varargin{1}))              % If input arguments are name-pair value
        obj.decide_which_state(varargin{:})                          % Ask proper function to create the right properties of the class
        
      elseif isnumeric(varargin{1}) && isnumeric(varargin{2})        % If input arguments are the first moments of the gaussian state
        R0 = varargin{1};
        V0 = varargin{2};
        
        assert(isvector(R0) && ismatrix(V0) && (length(R0) == length(V0)) && (size(V0, 1) == size(V0, 2)), "Unexpected first moments when creating gaussian state!") % Make sure they are a vector and a matrix with same length
        
        obj.R = R0(:);                 % Save mean quadratres   in a class properties (Parenthesis to ensure column vector)
        obj.V = V0;                    % Save covariance matrix in a class properties
        
      else
        error("Unexpected arguments when creating gaussian state!") % If input arguments do not make sense, call out the user
      end
      
    end
    
    function decide_which_state(obj, varargin)
      % If the user provided a name-pair argument to the constructor,
      % this function reads these arguments and creates the first moments of the gaussian state
      
      type_state = varargin{1};             % Name of expected type of gaussian state
      
      if strcmp(type_state, "vacuum")       % If it is a vacuum state
        obj.R = [0; 0];                     % Create its first moments
        obj.V = eye(2);
        return                              % End function
      end
      % Make sure there is an extra parameters that is a number
      assert(length(varargin)>1 && isnumeric(varargin{2}) && isscalar(varargin{2}), "Invalid or absent amplitude for non-vacuum gaussian state")
      
      if strcmp(type_state, "thermal")      % If it is a thermal state
        nbar = varargin{2};                 % Make sure its occuption number is a non-negative number
        assert(isreal(nbar) && (nbar>=0), "Imaginary or negative occupation number for thermal state")
        obj.R = [0; 0];
        obj.V = diag([2.0*nbar+1, 2.0*nbar+1]); % Create its first moments
        
      elseif strcmp(type_state, "coherent") % If it is a coherent state
        alpha = varargin{2};
        obj.R = [real(alpha); imag(alpha)];
        obj.V = eye(2);                     % Create its first moments
        
      elseif strcmp(type_state, "squeezed") % If it is a squeezed state
        r = varargin{2};                    % Make sure its squeezing parameter is a real number
        assert(isreal(r), "Unsupported imaginary amplitude for squeezed state")
        obj.R = [0; 0];
        obj.V = diag([exp(-2*r), exp(+2*r)]); % Create its first moments
        
      else
        error("Unrecognized gaussian state name, please check for typos or explicitelly pass its first moments as arguments")
      end
    end
    
    function rho = tensor_product(rho_A, rho_B)
      R0 = [rho_A.R; rho_B.R];
      V0 = blkdiag(rho_A.V, rho_B.V);
      
      rho = gaussian_state(R0, V0);
    end
    
    function rho_A = partial_trace(rho, indexes)
      N = length(rho.R)/2.0;
      N_A = length(rho.R) - 2.0*length(indexes); % Twice the number of modes in resulting state
      assert(N_A>=0, "Partial trace over more states than exists in gaussian state")
      
      modes = 1:N;
      entries = ~ismember(modes, indexes);
      modes = modes(entries);
      R0 = zeros(N_A, 1);
      for i=1:length(modes)
        j = modes(i);
        R0(2*i-1:2*i) = rho.R(2*j-1:2*j);
      end
      
      V0 = zeros(N_A, N_A);
      for i=1:length(modes)
        m = modes(i);
        for j=1:length(modes)
          n = modes(j);
          
          V0(2*i-1:2*i, (2*j-1:2*j)) = rho.V(2*m-1:2*m, 2*n-1:2*n);
        end
      end
      
      rho_A = gaussian_state(R0, V0);
    end
    
  end
end


%       R0 = -42*ones(N_A, 1);
%       for i=1:length(indexes)
%         j = indexes(i);
%         R0(2*i-1:2*i) = rho.R(2*j-1:2*j);
%       end

%elseif length(varargin)< 2 || ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
%  error("Invalid or absent amplitude for non-vacuum gaussian state")

% function obj = gaussian_state(varargin)
% default_coherent = 0;
% default_R0 = [0,0];
% default_V0 = eye(2);
% validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
% isvector(R0) && ismatrix(V0) && (length(R0) == length(V0))
% p = inputParser;
% validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%
% addOptional (p, 'R0', default_R0, @isvector);
% addOptional (p, 'V0', default_V0, @ismatrix);
% addParameter(p, 'coherent', default_coherent, @isstring);
% addParameter(p, 'shape',defaultShape,@(x) any(validatestring(x,expectedShapes)));
%
%
% parse(p,width,varargin{:});
% end
%
% end