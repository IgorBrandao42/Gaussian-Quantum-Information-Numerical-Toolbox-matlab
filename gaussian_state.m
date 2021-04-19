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
      
      if (isstring(varargin{1}) || ischar(varargin{1}))
        obj.decide_which_state(varargin{:})
        
      elseif isnumeric(varargin{1}) && isnumeric(varargin{2})
        R0 = varargin{1};
        V0 = varargin{2};
        if isvector(R0) && ismatrix(V0) && (length(R0) == length(V0))
          obj.R = R0;
          obj.V = V0;
        end
        
      else
        error("Unexpected first parameter when creating gaussian state!")
      end
      
    end
    
    function decide_which_state(obj, varargin)
      type_state = varargin{1};
      
      if strcmp(type_state, "vacuum")
        obj.R = [0, 0];
        obj.V = eye(2);
        return
      end
      
      assert(length(varargin)>1 && isnumeric(varargin{2}) && isscalar(varargin{2}), "Invalid or absent amplitude for non-vacuum gaussian state")
      
      if strcmp(type_state, "thermal")
        nbar = varargin{2};
        assert(isreal(nbar) && (nbar>=0), "Imaginary or negative occupation number for thermal state")
        obj.R = [0, 0];
        obj.V = diag([2.0*nbar+1, 2.0*nbar+1]);
        
      elseif strcmp(type_state, "coherent")
        alpha = varargin{2};
        obj.R = [real(alpha), imag(alpha)];
        obj.V = eye(2);
        
      elseif strcmp(type_state, "squeezed")
        r = varargin{2};
        assert(isreal(r), "Unsupported imaginary amplitude for squeezed state")
        
        obj.R = [0, 0];
        obj.V = diag([exp(-2*r), exp(+2*r)]);
        
      else
        error("Unrecognized gaussian state name, please check for typos or explicitelly pass its first moments as arguments")
      end
    end
    
  end
end


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