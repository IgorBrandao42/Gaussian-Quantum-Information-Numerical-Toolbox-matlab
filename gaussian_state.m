classdef gaussian_state < handle         % Class definning a nanoparticle
  properties
    R                                    % Mean quadratures
    V                                    % Covariance matrix
    
    Omega                                % Symplectic form matrix
    N_modes                              % Number of modes
  end
  
  methods
    % Constructor and its auxiliar functions
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
      %

      if nargin == 0                       % Default constructor (vacuum state)
        obj.R = [0; 0];                    % Save mean quadratres   in a class properties (Parenthesis to ensure column vector)
        obj.V = eye(2);                    % Save covariance matrix in a class properties
        obj.N_modes = 1;
      elseif (isstring(varargin{1}) || ischar(varargin{1}))              % If input arguments are name-pair value
        obj.decide_which_state(varargin{:})                          % Ask proper function to create the right properties of the class
        
      elseif isnumeric(varargin{1}) && isnumeric(varargin{2})        % If input arguments are the first moments of the gaussian state
        R0 = varargin{1};
        V0 = varargin{2};
        
        assert(isreal(R0) && isvector(R0) && ismatrix(V0) && (length(R0) == length(V0)) && (size(V0, 1) == size(V0, 2)), "Unexpected first moments when creating gaussian state!") % Make sure they are a vector and a matrix with same length
        
        obj.R = R0(:);                 % Save mean quadratres   in a class properties (Parenthesis to ensure column vector)
        obj.V = V0;                    % Save covariance matrix in a class properties
        obj.N_modes = length(R0)/2;
      else
        error("Unexpected arguments when creating gaussian state!") % If input arguments do not make sense, call out the user
      end
      
      omega = [[0, 1]; [-1, 0]];               % Auxiliar variable
      Omega = [];                              % Symplectic form matrix
      for i=1:obj.N_modes                      % Build the symplectic form
        Omega = blkdiag(Omega, omega);
      end
      obj.Omega = Omega;
      
      % Following line verify if the CM satisfy the uncertainty relation
      % obj.check_uncertainty_relation();
    end
    
    function V_check = check_uncertainty_relation(obj)
      % Check if the generated covariance matrix indeed satisfies the uncertainty principle
      
      V_check = obj.V + 1i*obj.Omega;
      assert(all(eig(V_check)>=0), "CM does not satisfy uncertainty relation!")
    end
    
    function decide_which_state(obj, varargin)
      % If the user provided a name-pair argument to the constructor,
      % this function reads these arguments and creates the first moments of the gaussian state
      
      obj.N_modes = 1;
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
        obj.R = [2*real(alpha); 2*imag(alpha)];
        obj.V = eye(2);                     % Create its first moments
        
      elseif strcmp(type_state, "squeezed") % If it is a squeezed state
        r = varargin{2};                    % Make sure its squeezing parameter is a real number
        assert(isreal(r), "Unsupported imaginary amplitude for squeezed state")
        obj.R = [0; 0];
        obj.V = diag([exp(-2*r), exp(+2*r)]); % Create its first moments
        
      else
        obj.N_modes = [];
        error("Unrecognized gaussian state name, please check for typos or explicitelly pass its first moments as arguments")
      end
    end
        
    % Construct another state, from this base gaussian_state
    function rho = tensor_product(rho_A, rho_array)
      % Given two gaussian states, calculates their tensor produtct, which is also a gaussian state
      
      % N_modes_array = [rho_array.N_modes];
      
      R_final = rho_A.R;                         % First moments of resulting state is the same of rho_A
      V_final = rho_A.V;                         % First block diagonal entry is the CM of rho_A
      
      for i=1:length(rho_array)
        R_final = cat(1,  R_final, rho_array(i).R);
        V_final = blkdiag(V_final, rho_array(i).V);
      end
      
      % Only valid for a tensor product with a single other state
      % R0 = [rho_A.R; rho_B.R];
      % V0 = blkdiag(rho_A.V, rho_B.V);
      
      rho = gaussian_state(R_final, V_final);
    end
    
    function rho_A = partial_trace(rho, indexes)
      % Partial trace over the modes in indexes of the complete gaussian state
      % The indexes are respective to the modes the user wants to trace out (as in the mathematical notation)
      
      N_A = length(rho.R) - 2.0*length(indexes); % Twice the number of modes in resulting state
      assert(N_A>=0, "Partial trace over more states than exists in gaussian state")
      
      % shouldnt there be an assert over max(indexes) < obj.N_modes ? -> you cant trace uot modes that do not exist
      
      modes = 1:rho.N_modes;
      entries = ~ismember(modes, indexes);
      modes = modes(entries);
      
      R0 = zeros(N_A, 1);
      V0 = zeros(N_A, N_A);
      
      for i=1:length(modes)
        m = modes(i);
        R0(2*i-1:2*i) = rho.R(2*m-1:2*m);
        
        for j=1:length(modes)
          n = modes(j);
          V0(2*i-1:2*i, (2*j-1:2*j)) = rho.V(2*m-1:2*m, 2*n-1:2*n);
        end
      end
      
      rho_A = gaussian_state(R0, V0);
    end
    
    function rho_A = only_modes(obj, indexes)
      % Partial trace over all modes except the ones in indexes of the complete gaussian state
      % The indexes are respective to the modes the user wants to retrieve
      
      N_A = length(indexes); % Twice the number of modes in resulting state
      assert(N_A>0 && N_A <= obj.N_modes, "Partial trace over more states than exists in gaussian state")
      
      R0 = zeros(N_A, 1);
      V0 = zeros(N_A, N_A);
      
      for i=1:length(indexes)
        m = indexes(i);
        R0(2*i-1:2*i) = obj.R(2*m-1:2*m);
        
        for j=1:length(indexes)
          n = indexes(j);
          V0(2*i-1:2*i, 2*j-1:2*j) = obj.V(2*m-1:2*m, 2*n-1:2*n);
        end
      end
      
      rho_A = gaussian_state(R0, V0);
    end
    
    % Properties of the gaussian state
    function p = purity(rho)
      % Purity of a gaussian state
      % Pure states have purity equal 1
      
      p = 1/prod(rho.symplectic_eigenvalues);
    end
    
    function [eta, V_sq, V_asq] = squeezing_degree(obj)
      % Degree of squeezing of the uncertainties of combinations of the quadratures
      % 
      % Calculates:
      % V_sq  - variance of the     squeezed quadrature
      % V_asq - variance of the antisqueezed quadrature
      % eta   - ratio of the variances above
      
      assert(obj.N_modes == 1, "At the moment, this program only calculates the wigner function for a single mode state")
      
      lambda = eig(obj.V);
      
      V_sq  = min(lambda);
      V_asq = max(lambda);
      
      eta = V_sq/V_asq;
    end
    
    function lambda = symplectic_eigenvalues(obj)
      % Calculates the sympletic eigenvalues of a covariance matrix V with symplectic form Omega
      % Finds the absolute values ofthe eigenvalues of i\Omega V and removes repeated entries
      
      H = 1i*obj.Omega*obj.V;                  % Auxiliar matrix
      lambda_0 = abs( eig(H) );                % Absolute value of the eigenvalues of the previous matrix
      
      lambda = zeros(obj.N_modes, 1);          % Variable to store the symplectic eigenvalues
      for i=1:obj.N_modes                      % Loop over the non-repeated entries of lambda_0
        lambda(i) = lambda_0(1);               % Get the first value on the repeated array
        lambda_0(1) = [];                      % Delete it
        
        [~, idx] = min( abs(lambda_0-lambda(i)) ); % Find the next closest value on the array (repeated entry)
        lambda_0(idx) = [];                    % Delete it too
      end
      
    end
    
    function Entropy = von_Neumann_Entropy(obj)
      % Calculation of the von Neumann entropy for a multipartite gaussian system
      %
      % MATHEMATICAL DESCRIPTION
      % Calculation of the von Neumann entropy for a multipartite gaussian system
      % given the covariance matrix that defines the state
      %
      % For any gaussian system, the von Neumann entropy (S_{N}) is a
      % function of the symplectic eigenvalues of the covariance matrix
      % (\nu_{k}):
      %
      % S_{N} = \sum_{k=1}^{N} g(\nu_k)
      %
      % where g(x) = [(x+1)/2]*log((x+1)/2) - [(x-1)/2]*log((x-1)/2) ,
      % and \nu_k are the sympletic eigenvalues of V,
      % i.e., modulus of the eigenvalues of i \Omega V
      
      nu = obj.symplectic_eigenvalues();  % Calculates the sympletic eigenvalues of a covariance matrix V
      
      nu(nu==1) = nu(nu==1) + 5e-16;             % 0*log(0) is NaN, but in the limit that x->0 : x*log(x) -> 0
                                          % MATLAB uses a 15 digits precision, I'm adding a noise at the limit of the numerical precision
                                          
      nu_plus  = (nu + 1)/2.0;            % Temporary variables
      nu_minus = (nu - 1)/2.0;
      g_nu = nu_plus.*log(nu_plus) - nu_minus.*log(nu_minus);
      
      Entropy = sum( g_nu );              % Calculate the entropy
    end
    
    function I = mutual_information(obj)
      % Mutual information for a multipartite gaussian system
      %
      % Calculates:
      % I     - mutual information  for the total system of the j-th covariance matrix
      % S_tot - von Neumann entropy for the total system of the j-th covariance matrix
      % S     - von Neumann entropy for the i-th mode    of the j-th covariance matrix
      
      S = zeros([obj.N_modes, 1]);               % Variable to store the entropy of each mode
      
      for j=1:obj.N_modes                        % Loop through each mode
        single_mode = obj.partial_trace(j);      % Get the covariance matrix for only the i-th mode
        S(j) = single_mode.von_Neumann_Entropy();% von Neumann Entropy for i-th mode of each covariance matrix
      end
      
      S_tot = obj.von_Neumann_Entropy();         % von Neumann Entropy for the total system of each covariance matrix
      
      I = sum(S) - S_tot;                        % Calculation of the mutual information
      
%       I = 0;                                     % Variable to store the mutual information
%       for j=1:obj.N_modes                        % Calculation of the mutual information
%         I = I + S(j);
%       end
%       I = I - S_tot;                             % Calculation of the mutual information
    end
    
    function W = wigner(obj, X, P)
      % Calculates the wigner function for the gaussian state
      % Currently only calculates for a single mode. 
      % More modes can be calculated, but will be difficult to visualize
      %
      % PARAMETERS
      %   X, P - 2D grid where the wigner function is to be evaluated
      
      assert(obj.N_modes == 1, "At the moment, this program only calculates the wigner function for a single mode state")
      
      % This squeezes the plotting in one direction and expands the other
      % No squeezing will apeear with this scaled axis !
%       if nargin == 1
%         x = obj.R(1) + 5*sqrt(obj.V(1,1))*linspace(-1, +1, 150);                 % Region to plot wigner function
%         p = obj.R(2) + 5*sqrt(obj.V(2,2))*linspace(-1, +1, 150);
%         [X, P] = meshgrid(x, p);
%       end
      
      N = obj.N_modes;                         % Number of modes
      W = zeros(length(X), length(P));         % Variable to store the calculated wigner function
      
      % I am not sure why there need an absolute value on W_den !
      % The theory does not need it, but numerical calculations were wrong only about the sign of the purity !
      
      one_over_purity = 1/obj.purity;
      
      for i=1:length(X)
        x = [X(i,:); P(i,:)];    
        for j=1:size(x, 2)
          dx = x(:, j) - obj.R;                % x_mean(:,i) is the i-th point in phase space
          
          W_num = exp( -(dx.')*(obj.V\dx)/2 ); % Numerator
          
          W_den = (2*pi)^N * one_over_purity;  % Denominator
          
          W(i, j) = W_num/W_den;               % Calculate the wigner function at every point on the grid
        end
      end
%       W = W.';
    end
    
    function F = fidelity(rho_1, rho_2)
      % Calculates the fidelity between the two arbitrary gaussian states
      % \rho_1 and \rho_2. See Phys. Rev. Lett. 115, 260501.
      %
      % The user should note that non-normalized quadratures are expected;
      % They are normalized to be in accordance with the notation of Phys. Rev. Lett. 115, 260501.
      
      assert(rho_1.N_modes == rho_2.N_modes, "Impossible to calculate the fidelity between gaussian states of diferent sizes!")
      
      u_1 = rho_1.R/sqrt(2.0);                   % Normalize the mean value of the quadratures
      u_2 = rho_2.R/sqrt(2.0);
      
      V_1 = rho_1.V/2.0;                         % Normalize the covariance matrices
      V_2 = rho_2.V/2.0;
      
      OMEGA = rho_1.Omega;
      
      delta_u = u_2 - u_1;                       % A bunch of auxiliar variables
      
      inv_V = inv(V_1 + V_2);
      
      V_aux = OMEGA.' * inv_V * ( OMEGA/4 + V_2*OMEGA*V_1 );
      
      identity = eye(2*rho_1.N_modes);
      
      F_tot_4 = det( 2*( sqrt(identity + (V_aux*OMEGA)^(-2)/4) + identity )*V_aux );
      
      F_0 = nthroot( real(F_tot_4) / det(V_1+V_2), 4); % We take only the real part of F_tot_4 as there can be a residual complex part from numerical calculations!
      
      F = F_0*exp( -delta_u.' * inv_V * delta_u  / 4); % Fidelity
      
    end
    
    function nbar = occupation_number(obj)
      % Occupation nuber for a single mode within the multipartite gaussian state
      
      Variances = diag(obj.V);              % From the current CM, take take the variances
      
      mean_x = obj.R(1:2:end);               % Odd  entries are the mean values of the position
      mean_p = obj.R(2:2:end);               % Even entries are the mean values of the momentum
      
      Var_x = Variances(1:2:end);              % Odd  entries are position variances
      Var_p = Variances(2:2:end);              % Even entries are momentum variances
      
      nbar = 0.25*( Var_x + mean_x.^2 + Var_p + mean_p.^2 ) - 0.5; % Calculate occupantion numbers at current time
      
    end
    
    % Entanglement
    function LN = logarithmic_negativity(obj, indexes)
      % Calculation of the logarithmic negativity for a bipartite system
      %
      % PARAMETERS:
      %    indexes - array with indices for the bipartition to consider 
      %    If the system is already bipartite, thos parameter is optional
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
      
      temp = obj.N_modes;
      if temp == 2                               % If the full system is only comprised of two modes
        V0 = obj.V;                              % Take its full covariance matrix
      elseif nargin > 1 && temp > 2
        assert(length(indexes) == 2, "Can only calculate the logarithmic negativity for a bipartition!");
        
        bipartition = obj.only_modes(indexes);   % Otherwise, get only the two mode specified by the user
        V0 = bipartition.V;                      % Take the full Covariance matrix of this subsystem
      end
      
      A = V0(1:2, 1:2);                       % Make use of its submatrices
      B = V0(3:4, 3:4);
      C = V0(1:2, 3:4);
      
      sigma = det(A) + det(B) - 2.0*det(C);  % Auxiliar variable
      
      ni = sigma/2.0 - sqrt( sigma^2 - 4.0*det(V0) )/2.0 ; % Square of the smallest of the symplectic eigenvalues of the partially transposed covariance matrix
      
      if ni < 0.0                            % Manually perform a maximum to save computational time (calculation of a sqrt can take too much time and deal with residual numeric imaginary parts)
        LN = 0.0;
      else
        ni = sqrt( real(ni) );               % Smallest of the symplectic eigenvalues of the partially transposed covariance matrix
        
        LN = max([0, -log(ni)]);             % Calculate the logarithmic negativity at each time
      end
    end
    
    % This is probabily wrong !
    function Duan = duan_criteria(obj, indexes) % I have to confirm the mathematics here!
      % Calculation of the LHS of the Duan criteria for a bipartite system
      %
      % PARAMETERS:
      %    indexes - array with indices for the bipartition to consider
      %    If the system is already bipartite, thos parameter is optional
      %
      % MATHEMATICAL DESCRIPTION
      % Calculation of the LHS of the Duan criteria for a bipartite system
      %
      % Given the quadratures (\hat{X}_j, \hat{P}_j) of the j-th subsystem of a bipartite system:
      %
      % \hat{X}_j =  (\hat{a]_j^\dagger + \hat{a}_j)  and
      % \hat{P}_j = i(\hat{a]_j^\dagger - \hat{a}_j)    ,
      %
      % the Duan criteria says that if a state is separable, then D > 1, where:
      %
      % D = 1/2 * { [\Delta(x_1 + x_2)]^2 + [\Delta(p_1 - p_2)]^2 } ,
      %
      % and ( \Delta(\hat{O}) )^2 = \langle \hat{O}^2 \rangle - \langle \hat{O} \rangle^2
      % is the variance of an operator \hat{0}.
      %
      % The Duan criteria is evaluated by use of the covariance matrix V respective to the 
      % bipartition: of modes indexes(1) and indexes(2) within the global covariance matrix of the state:
      % 
      % D =  1/2 * ( V(1,1) + V(2,2) + V(3,3) + V(4,4) + V(1,3) + V(3,1) - V(2,4) - V(4,2) )
      
      temp = obj.N_modes;
      if temp == 2                               % If the full system is only comprised of two modes
        V0 = obj.V;                              % Take its full covariance matrix
      elseif nargin > 1 && temp > 2
        assert(length(indexes) == 2, "Can only calculate the logarithmic negativity for a bipartition!");
        
        bipartition = obj.only_modes(indexes);   % Otherwise, get only the two mode specified by the user
        V0 = bipartition.V;                      % Take the full Covariance matrix of this subsystem
      end
      
    % Duan_0 = ( V0(1,1) + V0(2,2) + V0(3,3) + V0(4,4) )/4.0 + ( V0(1,3) - V0(2,4) )/2.0;
    % Duan_0 = Duan_0/2.0;
      
      Duan = ( V0(1,1) + V0(2,2) + V0(3,3) + V0(4,4) + V0(1,3) + V0(3,1) - V0(2,4) - V0(4,2) )/4;
    end
    
    % Gaussian unitaries (only applicable to a single mode states)
    function displace(obj, alpha)
      % Apply displacement operator on a single mode gaussian state
      % TO DO: generalize these operation to many modes!
    
      assert(     obj.N_modes   == 1, "Can only apply displacement operator on single mode state")
      d = [real(alpha); imag(alpha)];
      obj.R = obj.R + d;
      
      % If a displacement is attempted at a whole array of states, it is possible to apply a displacement in every entry
      % however, I cannot see why this would be the desired effect, I prefer to consider an error
      % assert(all([obj.N_modes]) == 1, "Can only apply displacement operator on single mode state")
      
    end
    
    function squeeze(obj, r)
      % Apply squeezing operator on a single mode gaussian state
      % TO DO: generalize these operation to many modes!
      
      assert(obj.N_modes == 1, "Error with input arguments when trying to apply displacement operator")
      S = diag([exp(-r), exp(+r)]);
      
      obj.R = S*obj.R;
      obj.V = S*obj.V*S;
    end
    
    function rotate(obj, theta)
      % Apply phase rotation on a single mode gaussian state
      % TO DO: generalize these operation to many modes!
      
      assert(obj.N_modes == 1, "Error with input arguments when trying to apply displacement operator")
      Rot = [[cos(theta), sin(theta)]; [-sin(theta), cos(theta)]];
      
      obj.R = Rot*obj.R;
      obj.V = Rot*obj.V*(Rot.');
    end
    
    % Gaussian unitaries (only applicable to a two mode states)
    function beam_splitter(obj, tau)
      % Apply a beam splitter transformation in a gaussian state
      % tau - transmissivity of the beam splitter
      
      assert(obj.N_modes==2, "Beam splitter transformation can only be applied for a two mode system");
      
      B = sqrt(tau)*eye(2);
      S = sqrt(1-tau)*eye(2);
      BS = [[B, S];[-S, B]];
      
      obj.R = BS*obj.R;
      obj.V = BS*obj.V*(BS.');
    end
    
    function two_mode_squeezing(obj, r)
      % Apply a two mode squeezing operator  in a gaussian state
      % r - squeezing parameter
      
      assert(obj.N_modes==2, "Two mode squeezing operator can only be applied for a two mode system");
      
      S0 = cosh(r)*eye(2);
      S1 = sinh(r)*diag([+  1,-1]);
      S2 = [[S0, S1];[S1, S0]];
      
      obj.R = S2*obj.R;
      obj.V = S2*obj.V*(S2.');
    end
    
  end % methods
end % classdef


% The following initialization of an array of gaussian states generates an unexpected behavior:
% coherent(1:3)  = gaussian_state("coherent", 2+1i);

% HOWEVER, the following works as intended
% coherent0  = gaussian_state("coherent", 2+1i);
% coherent = [coherent0, coherent0, coherent0];

% Both create an array of gaussian state, 
% but the first initialization seems to create a single object and store in every entry of the array
% So when a change is applied in one entry, the other ones are affected !