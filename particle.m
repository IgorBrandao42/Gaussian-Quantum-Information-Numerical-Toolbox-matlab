classdef particle < handle                                 % class definning a levitated nanoparticle
  properties
    omega                                                  % Natural frequency               [Hz] 
    g                                                      % Coupling strength with cavity   [Hz]
    gamma                                                  % Damping from the environment    [Hz]
    T_init                                                 % Initial temperature             [K]
    nbar_init                                              % Initial particle    occupation number 
    nbar_env                                               % Initial environment occupation number
    T_environment                                          % Temperature for its environment [K]
    
    nbar                                                   % Occupation number for each time
    Entropy                                                % Single mode entropy for each time
    is_entangled                                           % Boolean array telling at each time if it is entangled or not 
    T_approx                                               % Approximated/Effective temperature
    F_approx                                               % Fidelity between thermal state with approximated temperature and actual CM
    n_approx                                               % Approximated/Effective occupation number
  end
  
  methods
    function obj = particle(omega_0, g_0, gamma_0, T_init_0, T_environment_0)
      % Class simulating a nanoparticle interacting with an optical cavity field
      % The interaction is considered to be linear in the modes' positions quadratures
      % The particles is initially in a thermal state
      %
      % PARAMETERS:
      %   omega_0        - Natural frequencies
      %   g_0            - Coupling strength for the interaction with the cavity field
      %   gamma_0        - Mechanical damping 
      %   T_init_0       - Initial temperature
      %   T_environment_0 - Temperature for each heat bath
      %
      % Calculates the initial occupation number and occupation number of the associated heat bath
      
      obj.omega         = omega_0;
      obj.g             = g_0;
      obj.gamma         = gamma_0;
      obj.T_init        = T_init_0;
      obj.T_environment = T_environment_0;
      
      % Universal Constants
      k_B = 1.381e-23;                                     % Boltzmann's Constant [J/K]
      hbar = 1.055e-34;                                    % Planck's Constant    [J*s]
      
      obj.nbar_env  = 1/( exp(hbar*obj.omega./(k_B*obj.T_init)) - 1 ); 
      obj.nbar_init = 1/( exp(hbar*obj.omega./(k_B*obj.T_environment)) - 1 );
    end
    
    function approx_temperature(obj, V_0, k)
      % Internal function for 'simulation.fidelity_test'
      %
      % PARAMETERS:
      %   obj - class instance
      %   V_0 - single mode covariance matrix for this particle
      %   k   - index where to store approximated temperature
      
      if obj.is_entangled(k)                               % We can only approximate to a single moe thermal state, if it is not entangled!
        return
      end
      
      % Universal Constants
      k_B = 1.381e-23;                                     % Boltzmann's Constant [J/K]
      hbar = 1.055e-34;                                    % Planck's Constant    [J*s]
      
      options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000); % Options for the optimizer
      
      infidelity = @(Temperature) 1 - obj.fidelity_given_temperature(Temperature, V_0); % Define the infidelity function for the current particle's frequency
      
      T_min = 0;                                           % Minimum value the optimizer can return
      T_max = obj.T_init;                                  % Maximum value the optimizer can return
      
      obj.T_approx(k) = fminbnd(infidelity, T_min, T_max, options);         % Find the optimal temperature
        
      obj.F_approx(k) = obj.fidelity_given_temperature(obj.T_approx(k), V_0);  % Find the corresponding fidelity
      
      obj.n_approx(k) = 1/( exp( hbar*obj.omega/(k_B*obj.T_approx(k)) ) - 1 ); % Find the corresponding occupation number
    end
    
    function F = fidelity_given_temperature(obj, T, V_0)
      % Internal function for 'obj.approx_temperature'
      % 
      % Calculates the fidelity betwee the covariance matrix of a thermal state at temperature T
      % and the covariance matrix V_0
      % PARAMETERS:
      %   obj - class instance
      %   T   - temperature for thermal state
      %   V_0 - single mode covariance matrix to be compared
      
      k_B = 1.381e-23;                                 % J/K
      hbar = 1.055e-34;                                % J*s
      
      nbar_T = 1/( exp(hbar*obj.omega/(k_B*T)) - 1 );  % Occupation number for thermal at temperature T
      V_1 = diag([2*nbar_T + 1, 2*nbar_T + 1]);        % Covariance Matrix for a thermal state (zero mean)
      
      alpha = det(V_0 + V_1);                          % Auxiliar variables
      beta  = (det(V_0) - 1)*(det(V_1) - 1);           % Auxiliar variables
      
      F = 2.0/( sqrt(alpha+beta) - sqrt(beta) );       % Calculate the fidelity between the state with CM equal to V_0 and the thermal state at temperature T
      
    end
    
  end
end
  
  % A, D, V, V_ss, nbar, T_approx, F_approx, n_approx, J_particles, J_cavity, J_env, Log_Neg, Entropy, Entropy_single_mode
  

