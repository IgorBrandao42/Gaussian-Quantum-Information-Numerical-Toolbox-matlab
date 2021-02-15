classdef simulation < handle                               % Class simulating an experiment with N nanoparticles linearly interacting with a single optical cavity field
  properties                                               
    N_particles                                            % Number of particles
    Size_matrices                                          % Size of covariance, diffusion and drift matrices
    particles                                              % Array with particles (external class)
    cavity                                                 % Optical cavity       (external class)
    A                                                      % Drift matrix
    D                                                      % Diffusion Matrix
    V_0                                                    % Initial Covariance Matrix
    stable                                                 % Boolean telling if the system is stable or not
    
    mode_colors = [[0, 0, 0] ; [0, 0, 1]; [1, 0, 0]; [0, 0.5, 0]; [0.75, 0, 0.75]; [0.75, 0.75, 0]; [0.6350, 0.0780, 0.1840]]; % Array with colors to be consistently used throughout the plots
    
    t                                                      % Array with timestamps
    N_time                                                 % Length of time array
    Quadratures                                            % Semi classical time-evolved quadraures from Monte Carlos Langevin simulation
    V_cell                                                 % Cell with covariance matrix for each time
    V_ss                                                   % Steady-state Covariance Matrix calculated directly from Lyapunov equation
    J_particles                                            % Heat flux between every particle                        (expectation value)
    J_cavity                                               % Heat Flux between each  particle and the cavity         (expectation value)
    J_env                                                  % Heat Flux between each nanoparticle and the environment (expectation value, only valid when kappa > omega(j) )
    Log_Neg                                                % Logarithmic Negativity for each bipartition 
    Entropy_bipartitions                                   % von Neumann Entropy for each bipartition (see 'particle' class for single mode Entropy)
    Entropy_system                                         % von Neumann Entropy for the whole system
    I_tot                                                  % Mutual Information  for the whole system
    I_particles                                            % Mutual Information  for the particles
  end
  
  methods
    % TO DO: Generalize to more possible initial gaussian quantum states -> different inital conditions
    function obj = simulation(omega, g, gamma, T_init, T_environment, Delta, kappa)
      % Class simulating an experiment with N nanoparticles interacting with an optical cavity field.
      % The interaction is considered to be linear in the modes' positions quadratures
      % The cavity field is assumed to start in the vacuum state and the particles in a thermal state
      %
      % PARAMETERS:
      %   omega         - Array with Natural frequencies  of each nanoparticle
      %   g             - Array with Coupling strength   for each nanoparticle's interaction with the cavity field
      %   gamma         - Array with Mechanical damping  for each nanoparticle
      %   T_init        - Array with Initial temperature for each particle
      %   T_environment - Array with Temperature for each heat bath
      %   Delta         - Natural frequency of the optical mode (in the context of Coherent Scattering, this is the cavity-tweezer detunning)
      %   kappa         - Cavity linewidth
      %
      % The number of particles is inferred from the size of the input parameters (they must all be the same!)
      %
      % Calculates:
      % 'obj.A'      - the diffusion matrix of the system
      % 'obj.D'      - the drift matrix of the system
      % 'obj.V_0'    - the initial covariance matrix
      % 'obj.stable' - boolean telling if the system is stable or not
      
      number_parameters = [length(omega), length(g), length(gamma), length(T_init), length(T_environment)];
      if ~all( number_parameters == number_parameters(1))
        error("Wrong size of vectors with parameters! Impossible to determine number of particles, ending program...")
      end
      
      obj.cavity = optical_cavity(Delta, kappa);           % Optical cavity
      
      obj.N_particles = length(omega);                     % Number of particles
      obj.Size_matrices = 2*(obj.N_particles+1);           % Size of system and ccupation number for the environment (heat bath)
      
      obj.particles = {};                                  % Array with particles
      for j=1:obj.N_particles
        obj.particles{j} = particle(omega(j), g(j), gamma(j), T_init(j), T_environment(j));
      end
      
      % Universal Constants
      k_B = 1.381e-23;                                     % Boltzmann's Constant [J/K]
      hbar = 1.055e-34;                                    % Planck's Constant    [J*s]
      
      obj.A = [[-kappa/2 , Delta]; [-Delta   , -kappa/2]]; % Matrix defining the Langevin equation (without noises)
      for i=1:obj.N_particles
        A_particle = [[0, omega(i)]; [-omega(i), -gamma(i)]];
        obj.A = blkdiag(obj.A, A_particle);
        obj.A(2      , 2*i+1) = -2*g(i);
        obj.A(2*(i+1),   1    ) = -2*g(i);
      end
      
      nbar_environment = 1./( exp(hbar*omega./(k_B*T_environment)) - 1 );
      
      D_diag = zeros([obj.Size_matrices, 1]);              % Diagonal matrix with amplitude of noises autocorrelators
      D_diag(1:2) = [kappa; kappa];
      D_diag(4:2:end) = 2*gamma.*(2*nbar_environment + 1);
      obj.D = diag(D_diag);
      
      % Initial conditions for Lyapunov equation(vacuum state for cavity and thermal state for the partiles)
      nbar_0 = 1./( exp(hbar*omega./(k_B*T_init)) - 1 );   % Initial occupation numer for the particles
      obj.V_0 = [[1,0];[0,1]];                             % Covariance matrix for vacuum state
      for i=1:obj.N_particles                              % Covariance matrix for thermal satates
        obj.V_0 = blkdiag(obj.V_0, diag([2*nbar_0(i) + 1, 2*nbar_0(i) + 1]));
      end                                                  % The total state is separable and each mode has zero position and momentum average values, therefore blockdiagonal initial covariance matrix for the full system
      
      lambda = eig(obj.A);                                 % Calculate the eigenvalues of the diffusion matrix
      obj.stable = ~any( real(lambda) > 0 );               % Check if there is any with positive real part (unstability)
%       if ~obj.stable                                       % Warn the user about the stability
%         fprintf("The system is unstable! Matrix A has eigenvalues with positive real part !!\n\n");
%       else
%         fprintf("The system is stable!\n\n");
%       end
    end
    
    function run(obj, tspan, varargin)
      % Run every calculation on the simulation for the time stamps as input.
      % If additional parameters are passed, then it will only calculate what was specified.
      %
      % PARAMETERS:
      %   obj   - class instance
      %   tspan - Array with time stamps when the calculations should be done
      %   varargin (optional) - strings with the name of specific calculations that should be done
      %
      % Optional parameters are recommended to save computation time ! They are self-explanatory:
      % "occupation_number", "heat_flux", "entanglement", "entropy", "steady_state", "langevin", "fidelity_test"
      
      obj.t = tspan;                                       % Store the information about the time interval for the simulation 
      obj.N_time = length(tspan);                          % and its length on the simulation object
      
      if ~isempty(varargin)                                % If the user passed extra arguments
        what_to_calculate = string(varargin);
        obj.decide(what_to_calculate);                     % Only calculate what was asked for
      else
        obj.langevin();                                    % Calculate the semi-classical quadratures for each timestamp (high-temperature limit)
        
        obj.lyapunov();                                    % Calculate the CM for each timestamp (perform time integration of the Lyapunov equation)
        
        obj.occupation_number();                           % Calculate the occupation numbers for the cavity and particles
        
      % obj.heat_fluxes();                                 % Calculate the mean value of all heat fluxes
        
        obj.entanglement();                                % Calculate the logarithmic negativy for each bipartition of the system
        
        obj.entropy();                                     % Calculate the von Neumann entropy  for each bipartition of the system and every single mode von Neumann Entropy
        
        obj.steady_state();                                % Calculate steady-state CM and check if simulation achieved it
        
        obj.fidelity_test();                               % Approximate the effective temperature for each particle at each time
        
        obj.mutual_information();                          % Calculate the mutual information for the whole system
      end
      
    end
    
    function decide(obj, what_to_calculate)
      % If additional parameters were passed to function 'run', this function decides what the user asked for
      % It always perform 'lyapunov' method, without it the rest can not be calculated!
      %
      % PARAMETERS:
      %   obj   - class instance
      % 
      % This is an internal method and, in principle, does not need to be called by the user. Why are you reading this?!
      
      obj.lyapunov();                                    % Calculate the CM for each timestamp (perform time integration of the Lyapunov equation)
      
      if any(strcmp(what_to_calculate, "occupation_number"))
        obj.occupation_number();                         % Calculate the occupation numbers for the cavity and particles
      end
      
      if any(strcmp(what_to_calculate, "heat_flux"))
        obj.heat_fluxes();                               % Calculate the mean value of all heat fluxes
      end
      
      if any(strcmp(what_to_calculate, "entanglement"))
        obj.entanglement();                              % Calculate the logarithmic negativy  and von Neumann entropy for each bipartition of the system and every single mode von Neumann Entropy
      end
      
      if any(strcmp(what_to_calculate, "entropy"))
        obj.entropy();                                   % Calculate the logarithmic negativy  and von Neumann entropy for each bipartition of the system and every single mode von Neumann Entropy
      end
      
      if any(strcmp(what_to_calculate, "steady_state"))
        obj.steady_state();                              % Calculate steady-state CM and check if simulation achieved it
      end
      
      if any(strcmp(what_to_calculate, "langevin"))
        obj.langevin();                                  % Calculate the semi-classical quadratures for each timestamp (high-temperature limit)
      end
      
      if any(strcmp(what_to_calculate, "fidelity_test"))
        obj.fidelity_test();                             % Calculate approx. temperature for each particles at each timestamp through the Fidelity
      end
      
      if any(strcmp(what_to_calculate, "mutual_information"))
        obj.mutual_information();                        % Calculate the mutual information for the whole system at each timestamp
      end
      
    end
    
    % TO DO: Optimize use of classes
    function langevin(obj)
      % Solve the semi-classical Langevin equation for the expectation value of the quadrature operators
      % A Monte Carlos simulation to numericaly integrate the Langevin equations
      % The differential stochastic equations are solved through a Euler-Maruyama method
      %
      % The initial conditions for the particles follows the thermal state probability density in phase space
      % The initial conditions for the cavity    follows the vacuum  state probability density in phase space
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % Calculates 'obj.Quadratures', a matrix with the expectation values of the quadratures
      % obj.Quadratures(i,j) is the i-th quadrature expectation value at the j-th time
      %
      % Although correct, this method is poorly writen, optimizations to come!
      
      dt = obj.t(2) - obj.t(1);                            % Time step
      sq_dt = sqrt(dt);                                    % Square root of time step (for Wiener proccess in the stochastic integration)
      
      gamma    = zeros(obj.N_particles, 1);
      omega    = zeros(obj.N_particles, 1);
      T_init   = zeros(obj.N_particles, 1);
      nbar_env = zeros(obj.N_particles, 1);                % Store appropriately the environments occupation numbers
      for j=1:obj.N_particles
        nbar_env(j) = obj.particles{j}.nbar_init;
        T_init(j)   = obj.particles{j}.T_init;
        gamma(j)    = obj.particles{j}.gamma;
        omega(j)    = obj.particles{j}.omega;
      end
      
      N = zeros([obj.Size_matrices,1]);                    % Noise vector (fluctuation terms in Langevin equations)
      N(1) = obj.cavity.kappa; N(2) = N(1);
      N(4:2:end) = 2*gamma.*(2*nbar_env + 1);
      N = sqrt(N);
      
      % Parameters for the normal distribution of the initial conditions for Monte Carlos Simulations (vacuum state for cavity and thermal state for the partiles)
      % Universal Constants
      k_B = 1.381e-23;                                     % Boltzmann's Constant [J/K]
      hbar = 1.055e-34;                                    % Planck's Constant    [J*s]
      
      mu_cavity   = 0;                                     % Mean value for the position and momentum distribution for the vacuum  state
      mu_particle = 0;                                     % Mean value for the position and momentum distribution for the thermal state
      sigma_cavity = 1/sqrt(2);                            % Standard deviation of the position and momentum distribution for the vacuum state
      u = exp( - hbar*omega./(k_B*T_init) );               % Useful constant :)
      sigma_particle = sqrt( (1+u)./(1-u));                % Standard deviation of the position and momentum distribution for the thermal state (can be approximated by sqrt( 2*nbar_1 + 1 ) in the high temperature limit)
      
      obj.Quadratures = zeros([obj.Size_matrices, obj.N_time]); % Matrix to store each quadrature ensemble average at each time
      
    % rng(1)                                               % Default random number generator
      N_ensemble = 2e+2;                                   % Number of Monte Carlos iterations
    % milestone = 10; milestone_update = 10;               % Auxiliar variables to warn user that hthe computer didn't froze !
    % disp("Langevin simulation started...:")                 % Warn user that heavy calculation started
            
      for i=1:N_ensemble                                   % Loop on the random initial positions (% Monte Carlos simulation using Euler-Maruyama method in each iteration)
        
        X = zeros([obj.Size_matrices, obj.N_time]);        % For this iteration, this matrix stores each quadrature at each time (first and second dimensions, respectively)
        
        X(1,1) = normrnd(mu_cavity, sigma_cavity);         % Initial Cavity position quadrature (normal distribution for vacuum state)
        X(2,1) = normrnd(mu_cavity, sigma_cavity);         % Initial Cavity momentum quadrature (normal distribution for vacuum state)
        
        for j=1:obj.N_particles
          X(2*j+1, 1) = normrnd(mu_particle, sigma_particle(j)); % Initial Particle 1 position quadrature (normal distribution for thermal state)
          X(2*j+2, 1) = normrnd(mu_particle, sigma_particle(j)); % Initial Particle 1 momentum quadrature (normal distribution for thermal state)
        end
        
        for k = 1:obj.N_time-1                             % Euler-Maruyama method for stochastic integration
          X(:,k+1) = X(:,k) + (obj.A*X(:,k)*dt + sq_dt*N.*randn(size(N)));
        end % loop on time
        
      % milestone = update_user_on_loop_progress(i, N_ensemble, milestone, milestone_update); % Update the user on the calculation progress
        
        obj.Quadratures = obj.Quadratures + X;             % Add the new Monte Carlos iteration quadratures to the same matrix
        
      end % loop on each Monte Carlos realization
      
      obj.Quadratures = obj.Quadratures/N_ensemble;        % Divide the ensemble sum to obtain the average quadratures at each time
      
    % fprintf("Langevin simulation ended\n")             % Warn user that heavy calculation started
      
    end

    function lyapunov(obj)
      % Solve the lyapunov equation for the time evolved covariance matrix of the full system
      % 
      % Uses ode45 to numerically integrate, a fourth order Runge-Kutta method
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % Calculates 'obj.V_cell', a cell with the time evolved covariance matrix
      % obj.V_cell{j} is the covariance matrix at the j-th time
      
    % disp("Lyapunov simulation started...")               % Warn the user that heavy calculations started (their computer did not freeze!)
      
      V_0_vector = reshape(obj.V_0, [obj.Size_matrices^2, 1]); % Reshape the initial condition into a vector (expected input for ode45)
      
      [~, V_vector] = ode45(@(t,V) lyapunov_ode(t, V, obj.A, obj.D), obj.t, V_0_vector); % Solve Lyapunov eqaution through Fourth order Runge Kutta
      
      % Unpack the output of ode45 into a cell where each entry contains the information about the evolved CM at each time
      obj.V_cell = cell([size(V_vector,1),1]);                 % Initialize a cell to store all CMs for each time
      
      for i=1:size(V_vector,1)
        V_current_vector = V_vector(i, :);                 % Take the full Covariance matrix in vector form
        V_current = reshape(V_current_vector, [obj.Size_matrices, obj.Size_matrices]); % Reshape it into a proper matrix
        obj.V_cell{i} = V_current;                             % Store it in a cell
      end
      
   %  fprintf("Lyapunov simulation finished!\n\n")         % Warn user the heavy calculations ended
    end
    
    function occupation_number(obj)
      % Calculates the occupation number for each mode
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % Calculates the property 'obj.particles{i}.nbar(j)', a array with time evolved occupation number
      % for the i-th mode at the j-th time
      
      nbar = zeros([obj.N_particles+1, obj.N_time]);       % Initialize variable that will hold the occupation number for the cavity and particles
      
      for i=1:length(obj.V_cell)
        V = obj.V_cell{i};                                 % Take the full Covariance matrix at the current time
        Variances = diag(V);                               % From the current CM, take take the variances
        
        Var_x = Variances(1:2:end);                        % Odd  entries are position variances
        Var_p = Variances(2:2:end);                        % Even entries are momentum variances
        
        nbar(:, i) = 0.25*( Var_x + Var_p ) - 0.5;         % Calculate occupantion numbers at current time
      end
      
      obj.cavity.nbar = nbar(1, :);                        % Store the cavity     occupation number accordingly
      for j=1:obj.N_particles
        obj.particles{j}.nbar = nbar(j+1, :);              % Store the particles' occupation number accordingly
      end
      
    end
    
    function heat_fluxes(obj)
      % Calculate all the average heat fluxes present in the system, following arxiv:... (noy yet anounced)
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % Calculates the matrices 'obj.J_cavity', 'J_env' and 'J_particles':a cell with the time evolved covariance matrix
      % obj.V_cell{j} is the covariance matrix at the j-th time
      %
      % obj.J_cavity(i,j)     - Heat Flux from the cavity to the i-th particle at the j-th time
      % 
      % obj.J_env(i,j)        - Heat Flux from the environment to the i-th particle at the j-th time 
      % 
      % obj.J_particles(n,i,j)- Heat Flux between from the n-th particle to the i-th particle at the j-th time
      
      kappa = obj.cavity.kappa;
      Delta = obj.cavity.Delta;
      
      Omega = Delta/( (kappa/2)^2 + Delta^2);
      Gamma = Omega^2 * kappa/Delta;
      
      obj.J_particles = zeros([obj.N_particles, obj.N_particles, obj.N_time]); % Initialize variable that will hold the average heat flux between every particle
      obj.J_cavity    = zeros([obj.N_particles, obj.N_time]);                  % Initialize variable that will hold the average heat flux between each particle and the cavity
      obj.J_env       = zeros([obj.N_particles, obj.N_time]);                  % Initialize variable that will hold the average heat flux between each particle and the environment
      
      V_init = obj.V_cell{1};
      
      for i=1:length(obj.V_cell)                           % Loop through every timestamp
        V = obj.V_cell{i};                                 % Take the full Covariance matrix at the current time
        
        for n=1:obj.N_particles                            % Loop through every particle
          omega_n = obj.particles{n}.omega;
          g_n     = obj.particles{n}.g;
          gamma_n = obj.particles{n}.gamma;
          
          obj.J_cavity(n, i) = 2*omega_n*g_n*V(1, 2*n+2);   % Mean heat flux between cavity and n-th particle
          obj.J_env(n, i)    = omega_n*gamma_n*( V(2*n+2, 2*n+2) - V_init(2*n+2, 2*n+2) ); % Mean heat flux between environment and n-th particle
          
          for j=1:obj.N_particles                          % Loop through every particle
            omega_j = obj.particles{j}.omega;
            g_j     = obj.particles{j}.g;
                                                           % Mean heat flux between n-th particle and j-th particle
            obj.J_particles(n, j, i) = -4*pi*omega_j*g_j*g_n*( Gamma*omega_n*V(2*n+2,2*j+2) - Omega*V(2*n+1,2*j+2) ); 
          end % end loop over time
        end % end loop over particle n
      end % end loop over particle j
      
    end
    
    function entanglement(obj)
      % Calculates the logarithmic negativity for each bipartition
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % Calculates the property 'obj.Log_Neg', a array with time evolved logarithmic negativity
      % obj.Log_Neg(n, i, j) - Logarithmic Negativity between from the n-th particle and the i-th particle at the j-th time
      
      N_modes = obj.N_particles + 1;                       % Useful constant
      
      obj.Log_Neg = zeros([N_modes, N_modes, obj.N_time]); % Variable to store the Logarithmic Negativities
      
      for i=1:N_modes                                      % Loop through each mode
        for j=i+1:N_modes                                  % Loop through each non-repeated mode and skip diagonal
          
          V = bipartite_CM(obj.V_cell, i, j);          % From the full CM, get the sub matrix for the current bipartition
          
          obj.Log_Neg(i, j, :) = logarithmic_negativity2(V); % Logarithmic Negativity between modes i-1 and j-1
        end
      end                                                  % Mode 0 describes the cavity field and mode j>0 describe the j-th particle
      
    end
    
    function entropy(obj)
      % Calculates the von Neumann entropies for each bipartition and for each mode
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % Calculates the arrays:
      % obj.Entropy_bipartitions(n, i, j)      - Entropy in the bipartition with the n-th and i-th modes at the j-th time
      % obj.cavity.Entropy(j)       - Entropy for the cavity at the j-th time
      % obj.particles{n}.Entropy(j) - Entropy for the n-th particle at the j-th time
      
      N_modes = obj.N_particles + 1;                       % Useful constant
      
      obj.Entropy_bipartitions = zeros([N_modes, N_modes, obj.N_time]); % Variable to store the von Neumann Entropies
      Entropy_single_mode =    zeros([N_modes, obj.N_time]); % Variable to store the von Neumann Entropies
      
      for i=1:N_modes                                      % Loop through each mode
        V = single_mode_CM(obj.V_cell, i);                 % Get the CM for only the i-th mode
        
        Entropy_single_mode(i, :) = von_Neumann_Entropy(V);% von Neumann Entropy for i-th mode
        
        for j=i+1:N_modes                                  % Loop through each non-repeated mode
          V = bipartite_CM(obj.V_cell, i, j);          % From the full CM, get the sub matrix for the current bipartition
          
          obj.Entropy_bipartitions(i, j, :) = von_Neumann_Entropy(V); % von Neumann Entropy between modes i-1 and j-1
        end
      end                                                  % Mode 0 describes the cavity field and mode j>0 describe the j-th particle
      
      obj.cavity.Entropy = Entropy_single_mode(1, :).';      % Store the cavity entropy accordingly
      for j=1:obj.N_particles
        obj.particles{j}.Entropy = Entropy_single_mode(j+1, :).'; % Store the each particle entropy accordingly
      end
      
      obj.Entropy_system = von_Neumann_Entropy(obj.V_cell);% Calculate the entropy of the full system
      
    end
    
    function steady_state(obj)
      % Calculates the stedy state covariance matrix for the steady state
      % and checks if the final CM calculated is close enough to it (steady state achieved)
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % Calculates the matrix obj.V_ss with the steady state covariance matrix
      
      obj.V_ss = lyap(obj.A, obj.D);                       % Calculate steady-state CM
      
      V_ss_numerical = obj.V_cell{end};                    % Find the last CM found numerically
      V_diff = obj.V_ss - V_ss_numerical;
      relative_error = mean(abs(V_diff), 'all')./mean(abs(obj.V_ss), 'all');
      
      if max(relative_error) < 1e-2                        % Test if the system achieved the steady-state and warn user !
        fprintf("Steady state achieved!\n")
      else
        fprintf("Steady state has not be achieved!\n")
      end
      
    end
    
    % TO DO: Move from checking entanglement to checking the mutual information
    function fidelity_test(obj)
      % Calculates the best thermal state that approximates the covariance matrix of each particle at each time
      % Thermal state covariance matrices are compared to the single mode CM of each particle at each time
      % 
      % PARAMETERS:
      %   obj   - class instance
      %
      % Calculates:
      % obj.particles{i}.T_approx(j)     - Approximated temperature that approximates the covariance matrix of the i-th particle at the j-th time
      % obj.particles{i}.F_approx(j)     - Fidelity between corresponding thermal state for the i-th particle at the j-th time
      % obj.particles{i}.n_approx(j)     - Corresponding occupation number  for the i-th particle at the j-th time
      % obj.particles{i}.is_entangled(j) - Boolean telling if the for the i-th particle is entangled at the j-th time
      % The last varaible is DEPRECATED and will be removed
      
      if isempty(obj.Log_Neg)
        disp("No previous logarithmic negativity calculated, verifying if any nanoparticle became entangled...")
        obj.entanglement();
      end
      
      for j=1:obj.N_particles                                 % Loop through every particle
        obj.particles{j}.T_approx =     zeros(obj.N_time, 1); % Initialize everything to spurios values
        obj.particles{j}.F_approx =     zeros(obj.N_time, 1); % only populate it if the particle is disentangled
        obj.particles{j}.n_approx =     zeros(obj.N_time, 1);
        obj.particles{j}.is_entangled = false(obj.N_time, 1); % Tell if the j-th particle is entangled at each time
      end
      
      for k=1:obj.N_time
        V_full = obj.V_cell{k};                            % Take the full Covariance matrix at the current time
        currently_entangled_matrix = obj.Log_Neg(:,:,k)>0; % Get matrix telling if i-th particle is entangled with j-th particle
        
        for j=1:obj.N_particles                            % Loop through every particle
          entalged_b4  = currently_entangled_matrix(j+1, :); % Array telling if current particle is entangle with some other previous   mode
          entalged_next= currently_entangled_matrix(:, j+1); % Array telling if current particle is entangle with some other subsequent mode
          obj.particles{j}.is_entangled(k) = any(entalged_b4) | any(entalged_next);    % j+1 , because we have to skip the cavity !
          
          V_particle = V_full(2*j+1:2*j+2, 2*j+1:2*j+2);   % Single mode covariance matrix for the current particle
          
          obj.particles{j}.approx_temperature(V_particle, k); % Approximate its effective temperature through the Fidelity
          
        end % end of time loop
      end % end of particles loop
      
    end
    
    function mutual_information(obj)
      % Calculates the mutual information for the complete system
      % 
      % PARAMETERS:
      %   obj   - class instance
      % 
      % Calculates:
      % obj.I_tot(j) - mutual information for the complete system at the j-th time
      
      if isempty(obj.particles{1}.Entropy) || isempty(obj.Entropy_bipartitions)
        obj.entropy();                                     % If no entropy was previously calculated, calculate it!
      end
      
      obj.I_tot = zeros( size(obj.cavity.Entropy) );                   % We start with the entropy of the cavity
      
      for j=1:obj.N_particles
        obj.I_tot = obj.I_tot + obj.particles{j}.Entropy; % Then we sum with the entropy of each nanoparticle
      end
      S_each_particle_sum = obj.I_tot;                    % Save the sum of the particles' entropy for latter
      
      obj.I_tot = obj.I_tot + obj.cavity.Entropy;         % Sum up the final mode entropy
      
      obj.I_tot = obj.I_tot - obj.Entropy_system;         % And finally subtract the entropy of the full system
      
      V_particles = obj.CM_for_the_particles();           % Find the covariance matrix just for the particles
      S_particles = von_Neumann_Entropy(V_particles);     % Find the entropy for all particles together
      
      obj.I_particles = S_each_particle_sum - S_particles;%Mutual information for the particles
    end
    
    function V = CM_for_the_particles(obj)
        % Finds the covariance submatrices only for the particles
        %
        % PARAMETERS:
        %   obj   - class instance
        
        
        V = cell( size(obj.V_cell) );                % Cell of the same size, but each entry is a single mode CM
        V_aux = obj.V_cell;
        
        parfor i=1:length(obj.V_cell)
          V_temp = V_aux{i};                    % Take the full Covariance matrix at the i-th entry
          
          V{i} = V_temp(3:end, 3:end);           % Only look for the submatrix corresponding to the desired mode
        end
        
    end
    
    function plot_entanglement(obj)
      if isempty(obj.Log_Neg)
        disp("No logarithmic negativity, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Logarithmic Negativity for each bipartition"); % Open a figure to contain the subplots
      
      title_begining = "Particle " + string(1:obj.N_particles); % Create the begining of the title for each mode
      title_begining = ["Cavity", title_begining];
      
      N_modes = obj.N_particles + 1;                       % Useful constant
      
      N_bipartitions = N_modes*(N_modes-1)/2;               % Number of entries above the diagonal
      nx = divisors(N_bipartitions);                        % Find out how to subdivide the plotting figure
      if length(nx)>1                                       % into many subplots
        nx = nx(2);
      end
      ny = N_bipartitions/nx;
      idx = 1;
      
      for i=1:N_modes                                      % Loop through each mode
        for j=i+1:N_modes                                  % Loop through each non-repeated mode
          subplot(nx, ny, idx)                  % Create subplots on the same figure
          idx = idx + 1;
          
          plot(obj.t, reshape(obj.Log_Neg(i, j, :), [obj.N_time, 1]), 'b')
          ylabel("Log. Neg.");
          
          title(title_begining(i) + " and particle " + (j-1));          % Title for current subplot
          
          xlim([obj.t(1), obj.t(end)]);                                 % Adjust x-axis and its label
        % xticks((0:N_oscillations/5:N_oscillations)/OMEGA);
        % x_label = string(0:N_oscillations/5:N_oscillations);
        % xticklabels(x_label)
          xlabel('t [s]')
          
          set(gca, 'Fontsize', 18)                                       % Increase the font for readability
        end
      end
    end
    
    function plot_entanglement_and_entropy(obj)
      if isempty(obj.Log_Neg) || isempty(obj.Entropy_bipartitions)
        disp("No logarithmic negativity or entropy calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Logarithmic Negativity and von Neumann Entropy for each bipartition"); % Open a figure to contain the subplots
      
      title_begining = "Particle " + string(1:obj.N_particles); % Create the begining of the title for each mode
      title_begining = ["Cavity", title_begining];
      
      N_modes = obj.N_particles + 1;                       % Useful constant
      
      N_bipartitions = N_modes*(N_modes-1)/2;               % Number of entries above the diagonal
      nx = divisors(N_bipartitions);                        % Find out how to subdivide the plotting figure
      if length(nx)>1                                       % into many subplots
        nx = nx(2);
      end
      ny = N_bipartitions/nx;
      idx = 1;
      
      for i=1:N_modes                                      % Loop through each mode
        for j=i+1:N_modes                                  % Loop through each non-repeated mode
          subplot(nx, ny, idx)                  % Create subplots on the same figure
          idx = idx + 1;
          
          yyaxis left                                                   % Left axis of the plot is entanglement
          plot(obj.t, reshape(obj.Log_Neg(i, j, :), [obj.N_time, 1]), 'b')
          ylabel("Log. Neg.");
          
          yyaxis right                                                  % Right axis of the plot is entropy
          plot(obj.t, reshape(obj.Entropy_bipartitions(i, j, :), [obj.N_time, 1]), 'r')
          ylabel("Entropy");
          
          title(title_begining(i) + " and particle " + (j-1));          % Title for current subplot
          
          xlim([obj.t(1), obj.t(end)]);                                 % Adjust x-axis and its label
        % xticks((0:N_oscillations/5:N_oscillations)/OMEGA);
        % x_label = string(0:N_oscillations/5:N_oscillations);
        % xticklabels(x_label)
          xlabel('t [s]')
          
          set(gca, 'Fontsize', 18)                                       % Increase the font for readability
        end
      end
    end
    
    function plot_single_mode_entropy(obj)
      if isempty(obj.cavity.Entropy)
        disp("No single mode entropy calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Single mode von Neumann Entropy");   % Open a figure to contain the subplots
      hold on
      
      title_begining = "Particle " + string(1:obj.N_particles); % Create the begining of the title for each mode
      
      plot(obj.t, obj.cavity.Entropy, 'Linewidth', 1.5, 'DisplayName', "Cavity", 'Color', obj.mode_colors(1,:))
      for k=1:obj.N_particles                              % Calculate and plot entanglement and entropy
        plot(obj.t, obj.particles{k}.Entropy, 'Linewidth', 1.5, 'DisplayName', title_begining(k), 'Color', obj.mode_colors(k+1,:))
      end
      
      ylabel("Entropy");
      title("Single mode von Neumann Entropy");             % Title for current subplot
      
      xlim([obj.t(1), obj.t(end)]);                         % Adjust x-axis and its label
      xlabel('t [s]')
      
      legend('Location', 'best')                           % Show the legend
      set(gca, 'Fontsize', 18)                             % Increase the font for readability
      set(gca, 'TickDir', 'out')                           % Set thicks outward
    end
    
    function plot_system_entropy(obj)
      if isempty(obj.Entropy_system)
        disp("The entropy of the system was not calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Complete system's von Neumann Entropy");      % Open a figure to contain the subplots
      hold on
      
      title("Complete system's von Neumann Entropy");     % Title for the plot
      
      plot(obj.t, real(obj.Entropy_system), 'k', 'Linewidth', 1.5)
      
      ylabel("System Entropy");
      
      xlim([obj.t(1), obj.t(end)]);                        % Adjust x-axis and its label
      xlabel('t [s]')
      
      set(gca, 'Fontsize', 18)                             % Increase the font for readability
      set(gca, 'TickDir', 'out')                           % Set thicks outward
    end
    
    function plot_mutual_information(obj)
      if isempty(obj.I_tot)
        disp("The mutual information of the system was not calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Mutual information of the system");      % Open a figure to contain the subplots
      hold on
      
      title("Mutual information of the system");     % Title for the plot
      
      plot(obj.t, real(obj.I_tot), 'k', 'Linewidth', 1.5)
      
      ylabel("Mutual information");
      
      xlim([obj.t(1), obj.t(end)]);                        % Adjust x-axis and its label
      xlabel('t [s]')
      
      set(gca, 'Fontsize', 18)                             % Increase the font for readability
      set(gca, 'TickDir', 'out')                           % Set thicks outward
      
      if isempty(obj.I_particles)
        disp("The mutual information of the particles was not calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Mutual information of the particles");      % Open a figure to contain the subplots
      hold on
      
      title("Mutual information of the particles");     % Title for the plot
      
      plot(obj.t, real(obj.I_particles), 'k', 'Linewidth', 1.5)
      
      ylabel("Mutual information");
      
      xlim([obj.t(1), obj.t(end)]);                        % Adjust x-axis and its label
      xlabel('t [s]')
      
      set(gca, 'Fontsize', 18)                             % Increase the font for readability
      set(gca, 'TickDir', 'out')                           % Set thicks outward
    end
    
    function plot_occupation_number(obj)
      if isempty(obj.cavity.nbar)
        disp("No occupation number calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Occupation number for each mode"); % Open a figure to contain the subplots
      hold on
      title("Occupation number for each mode")
      
      title_begining = "Particle " + string(1:obj.N_particles); % Create the begining of the title for each mode
      
      plot(obj.t, obj.cavity.nbar, 'Linewidth', 1.5, 'DisplayName', "Cavity", 'Color', obj.mode_colors(1,:))
      for k=1:obj.N_particles
        plot(obj.t, obj.particles{k}.nbar, 'Linewidth', 1.5, 'DisplayName', title_begining(k), 'Color', obj.mode_colors(k+1,:))
      end
      
      ylabel('$\bar{n}$','Interpreter','Latex','rotation',0,'VerticalAlignment','middle');
      xlabel('t [s]');
      xlim([obj.t(1), obj.t(end)]);
      
    % xticks((0:N_oscillations/5:N_oscillations)/OMEGA);
    % x_label = string(0:N_oscillations/5:N_oscillations);
    % xticklabels(x_label)
      
      legend
      set(gca, 'Fontsize', 18)
      set(gca, 'TickDir', 'out')
      
      set(gca, 'YScale', 'log')
    end
    
    function plot_heat_fluxes(obj)
      if isempty(obj.J_cavity)
        disp("No average heat flux calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Heat Fluxes between cavity and nanoparticles"); % Open a figure to contain the subplots
      hold on
      title("Heat Fluxes between cavity and nanoparticles")
      
      for j=1:obj.N_particles
        plot(obj.t, obj.J_cavity(j,:), 'Linewidth', 2, 'DisplayName', "cav, " + j, 'Color', obj.mode_colors(j+1,:));
      end
      
      xlim([obj.t(1), obj.t(end)]);                        % Adjust x-axis and its label
      xlabel('t [s]')
      legend('interpreter', 'latex', 'Location', 'best')
      set(gca, 'Fontsize', 18)                             % Increase the font for readability
      ylabel("J_{cavity, j} [a.u.]")
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      figure('Name', "Heat Fluxes between cavity and nanoparticles for"); % Open a figure to contain the subplots
      hold on
      title("Heat Fluxes between environment and nanoparticles")
      
      for j=1:obj.N_particles
        plot(obj.t, obj.J_env(j,:), 'Linewidth', 2, 'DisplayName', "env, " + j, 'Color', obj.mode_colors(j+1,:));
      end
      
      xlim([obj.t(1), obj.t(end)]);                        % Adjust x-axis and its label
      xlabel('t [1s]')
      legend('interpreter', 'latex', 'Location', 'best')
      set(gca, 'Fontsize', 18)                             % Increase the font for readability
      ylabel("J_{env, j} [a.u.]") 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      figure('Name', "Heat Fluxes between nanoparticles"); % Open a figure to contain the subplots
      hold on
      title("Heat Fluxes between particles")
      
      idx = 2;
      for n=1:obj.N_particles
        for j=n+1:obj.N_particles
          plot(obj.t, reshape(obj.J_particles(n,j,:), obj.N_time, 1), 'Linewidth', 2, 'DisplayName', n + "," + j);
          idx = idx + 1;
        end
      end
      
      xlim([obj.t(1), obj.t(end)]);                             % Adjust x-axis and its label
      xlabel('t [1s]')
      legend('interpreter', 'latex', 'Location', 'best')
      set(gca, 'Fontsize', 18)                          % Increase the font for readability
      ylabel("J_{particles} [a.u.]") % [(2\omega_1)]") % /(2*omega(1))
      
    end
    
    function plot_phase_space(obj)
      if isempty(obj.Quadratures)
        disp("No mean quadratures calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Quadratures Phase Space");           % Open a figure to contain the subplots
      
      x_label = "x_" + string(1:obj.N_particles);          % Create the labels for the x axis
      x_label = ["q", x_label];
      
      y_label = "p_" + string(1:obj.N_particles);          % Create the labels for the y axis
      y_label = ["p", y_label];
      
      title_name = "Particle " + string(1:obj.N_particles);% Create the title (each mode name)
      title_name = ["Cavity", title_name];
      
      for i=0:obj.N_particles                              % For each mode
        x = obj.Quadratures(2*i+1, :);                     % Position quadrature for the i-th mode
        p = obj.Quadratures(2*i+2, :);                     % Momentum quadrature for the i-th mode
        
        dq = x(2:end)-x(1:end-1);                          % Increments in the position for each time step to be used in the arrow plot
        dp = p(2:end)-p(1:end-1);                          % Increments in the momentum for each time step to be used in the arrow plot
        
        quiver_scale = 0;                                  % Tell MATLAB to not scale the arrows
        
        subplot(1, obj.N_particles+1, i+1);                % Specify why subplot to use
        quiver(x(1:end-1), p(1:end-1), dq, dp, quiver_scale, 'Color', obj.mode_colors(i+1,:)) % Arrow plot the phase space trajectory
        hold on
        plot(x(1), p(1), 'x', 'Color', obj.mode_colors(i+1,:), 'MarkerSize', 12)
        
        title(title_name(i+1))                             % Add a title
        xlabel(x_label(i+1));                              % Label the x axis
        ylabel(y_label(i+1),'rotation',0,'VerticalAlignment','middle'); % Label the p axis
        
        set(gca, 'Fontsize', 18)                           % Set text to be on proper font size
        
      end
      
    end
    
    function plot_fidelity_approx(obj)
      if isempty(obj.particles{1}.T_approx)
        disp("No temperature approximation calulated, please calculate it before plotting it!")
        return
      end
      
      figure('Name', "Temperature approximation");         % Open a figure to contain the subplots
      hold on
      title_name = "Particle " + string(1:obj.N_particles);% Create the title (each mode name)
      
      for j=1:obj.N_particles                              % For each particle
        subplot(obj.N_particles, 1, j);                    % Specify why subplot to use
        
        yyaxis left                                        % Left axis of the plot is approx. temperature
        plot(obj.t, obj.particles{j}.T_approx, 'b-', 'Linewidth', 2); %, 'Color', obj.mode_colors(j+1,:)) % Plot approximate temperature
        ylabel("T_{approx} [K]", 'rotation', 90,'VerticalAlignment','middle'); % Label the left y axis
        xlim([obj.t(1), obj.t(end)]);                      % Adjust x-axis and its label
      % set(gca, 'YScale', 'log')
        
        yyaxis right                                       % Right axis of the plot is fidelity of temperature
        plot(obj.t, obj.particles{j}.F_approx, 'r-', 'Linewidth', 2)        % Plot the fidelity of the approximatioin
        ylabel("F_{approx}",'rotation',90,'VerticalAlignment','middle'); % Label the left y axis
        ylim([0.0, 1.0])
        
        title(title_name(j))                               % Add a title
        xlabel("t [s]");                                   % Label the x axis
        xlim([obj.t(1), obj.t(end)]);                      % Adjust x-axis and its label
        set(gca, 'Fontsize', 18)                           % Set text to be on proper font size
      end
      
    end
    
    function plot(obj)
      obj.plot_entanglement_and_entropy()
      
    % obj.plot_entanglement();
      
      obj.plot_single_mode_entropy()
      
    % obj.plot_heat_fluxes()
      
      obj.plot_occupation_number()
      
      obj.plot_phase_space();
      
      obj.plot_fidelity_approx();
      
      obj.plot_system_entropy();
      
      obj.plot_mutual_information();
    end
    
  end % end methods
end % end class




