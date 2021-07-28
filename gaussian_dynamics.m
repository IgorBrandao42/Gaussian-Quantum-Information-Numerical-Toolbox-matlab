classdef gaussian_dynamics < handle              % Class simulating the time evolution of a gaussian state following a set of Langevin and Lyapunov equations for its first moments dynamics
  properties
    A                                            % Drift matrix (can be a function handle to have a time dependency!)
    D                                            % Diffusion Matrix 
    N                                            % Mean values of the noises
    initial_state                                % Initial state of the global system
    t                                            % Array with timestamps for the time evolution
    
    N_time                                       % Length of time array
    Size_matrices                                % Size of covariance, diffusion and drift matrices
    
    R_semi_classical                             % Array with semi-classical mean quadratures (Semi-classical time evolution using Monte Carlos method)
    R                                            % Array with mean quadratures  for each time
    V                                            % Cell  with covariance matrix for each time
    state                                        % Gaussian state               for each time
    
    is_stable                                    % Boolean telling if the system is stable or not
    steady_state_internal                        % Steady state
  end
  
  methods
    function obj = time_evolution(A_0, D_0, N_0, initial_state_0)
      % Class constructor for simulating the time evolution of the global system
      % Open/closed quantum dynamics dictated by Langevin and Lyapunov equations
      %
      % Langevin: \dot{R} = A*X + N           : time evolution of the mean quadratures
      %
      % Lyapunov: \dot{V} = A*V + V*A^T + D   : time evolution of the covariance matrix
      %
      % PARAMETERS:
      %    A_0           - Drift matrix     (numerical matrix or funciton handle for a matrix with time dependency
      %    D_0           - Diffusion Matrix (auto correlation of the noises, assumed to be delta-correlated in time)
      %    N_0           - Mean values of the noises
      %    initial_state - Cavity linewidth
      %
      % CALCULATES:
      %    obj           - instance of a time_evolution class
      %    obj.is_stable - boolean telling if the system is stable or not
      
      obj.A = A_0;                               % Drift matrix
      obj.D = D_0;                               % Diffusion Matrix
      obj.N = N_0;                               % Mean values of the noises
      
      obj.initial_state = initial_state_0;       % Initial state of the global system
      
      obj.Size_matrices = length(obj.D);         % Size of system and ccupation number for the environment (heat bath)
      
    % assert(2*initial_state_0.N_modes == obj.Size_matrices); % Check if the initial state and diffusion/drift matrices have appropriate sizes !
      
      if ~isa(obj.A,'function_handle')
        lambda = eig(obj.A);                     % Calculate the eigenvalues of the drift matrix
        obj.is_stable = ~any( real(lambda) > 0 );% Check if there is any with positive real part (unstability)
      end
    end
    
    function result = run(obj, t_span)
      % Run every time evolution available at the input timestamps.
      %
      % PARAMETERS:
      %   obj   - class instance
      %   tspan - Array with time stamps when the calculations should be done
      %
      % CALCULATES:
      %   result = array with time evolved gaussian states for each timestamp of the input argument t_span
      %   each entry of the array is a gaussian_state class instance
      
      obj.langevin(t_span);                      % Calculate the mean quadratures for each timestamp
      
      obj.lyapunov(t_span);                      % Calculate the CM for each timestamp (perform time integration of the Lyapunov equation)
      
      obj.build_states();                        % Combine the time evolutions calculated above into an array of gaussian states
      
      result = obj.state;
      
    % obj.langevin_semi_classical();             % Calculate the semi-classical quadratures for each timestamp
      
    % obj.steady_state();                        % Calculate the gaussian steady state of the system
      
    end
    
    function langevin(obj, t_span)
      % Solve the Langevin equation for the time evolved mean quadratures of the full system
      %
      % Uses ode45 to numerically integrate the average Langevin equations (a fourth order Runge-Kutta method)
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % CALCULATES:
      %   obj.R - a cell with the time evolved mean quadratures where
      %   obj.R(i,j) is the i-th mean quadrature at the j-th timestamp
      
      obj.t = t_span;                                 % Timestamps for the simulation
      obj.N_time = length(t_span);                    % Number of timestamps
      
      if isa(obj.A,'function_handle')                 % I have to check if there is a time_dependency on the odes :(
        langevin_ode = @(t,R) obj.A(t)*R + obj.N;     % Function handle that defines the Langevin equation (returns the derivative)
      else
        langevin_ode = @(t,R) obj.A*R + obj.N;        % Function handle that defines the Langevin equation (returns the derivative)
      end
      
      [~, R_temp] = ode45(langevin_ode, obj.t, obj.initial_state.R);   % Solve Langevin eqaution through Fourth order Runge Kutta
      % Each row in R corresponds to the solution at the value returned in the corresponding row of obj.t
      
      obj.R = R_temp.';                               % Initialize a cell to store all vetor of quadratures for each time
      
      %  fprintf("Langevin simulation finished!\n\n") % Warn user the heavy calculations ended
    end
    
    function lyapunov(obj, t_span)
      % Solve the lyapunov equation for the time evolved covariance matrix of the full system
      %
      % Uses ode45 to numerically integrate, a fourth order Runge-Kutta method
      %
      % PARAMETERS:
      %   obj   - class instance
      %
      % CALCULATES:
      %  'obj.V' - a cell with the time evolved covariance matrix where
      %   obj.V{j} is the covariance matrix at the j-th timestamp
      
      % disp("Lyapunov simulation started...")             % Warn the user that heavy calculations started (their computer did not freeze!)
      
      obj.t = t_span;                                      % Timestamps for the simulation
      obj.N_time = length(t_span);                         % Number of timestamps
      
      V_0_vector = reshape(obj.initial_state.V, [obj.Size_matrices^2, 1]); % Reshape the initial condition into a vector (expected input for ode45)
      
      if isa(obj.A,'function_handle')                      % I have to check if there is a time_dependency on the odes :(
        ode = @(t,V) lyapunov_ode(t, V, obj.A(t), obj.D);  % Function handle that defines the Lyapunov equation (returns the derivative)
      else
        ode = @(t,V) lyapunov_ode(t, V, obj.A, obj.D);     % Function handle that defines the Lyapunov equation (returns the derivative)
      end
      
      [~, V_vector] = ode45(ode, obj.t, V_0_vector);       % Solve Lyapunov equation through Fourth order Runge Kutta
      
      % Unpack the output of ode45 into a cell where each entry contains the information about the evolved CM at each time
      obj.V = cell([size(V_vector,1),1]);                  % Initialize a cell to store all CMs for each time
      
      for i=1:size(V_vector,1)
        V_current_vector = V_vector(i, :);                 % Take the full Covariance matrix in vector form
        V_current = reshape(V_current_vector, [obj.Size_matrices, obj.Size_matrices]); % Reshape it into a proper matrix
        obj.V{i} = V_current;                              % Store it in a cell
      end
      
      %  fprintf("Lyapunov simulation finished!\n\n")      % Warn user the heavy calculations ended
    end
    
    function build_states(obj)
      % Builds the gaussian state at each time from their mean values and covariance matrices
      % This funciton is completely unnecessary, but it makes the code more readable :)
      %
      % CALCULATES:
      %   obj.state - array with time evolved gaussian states for each timestamp of the input argument t_span
      %   each entry of the array is a gaussian_state class instance
      
      assert(~isempty(obj.R) && ~isempty(obj.V), "No mean quadratures or covariance matrices, can not build time evolved states!")
      
      temp = repmat(gaussian_state(), obj.N_time, 1);       % Initialize an "empty" array of gaussian states
      
      RR = obj.R;
      VV = obj.V;
      
      for i=1:obj.N_time
        temp(i) = gaussian_state(RR(:, i), VV{i});
      end
      obj.state = temp;
    end
    
    function ss = steady_state(obj, A_0, A_c, A_s, omega)
      % Calculates the steady state for the system
      %
      % PARAMETERS:
      %   obj   - class instance
      % 
      %   The next parameters are only necessary if the drift matrix has a time dependency (and it is periodic)
      %   A_0, A_c, A_s - components of the Floquet decomposition of the drift matrix
      %   omega - Frequency of the drift matrix
      % 
      % CALCULATES:
      %   obj.steady_state_internal with the steady state (gaussian_state)
      %   ss - gaussian_state with steady state of the system
      
      if ~isa(obj.A, 'numeric')                  % If the Langevin and Lyapunov eqs. have a time dependency, move to the Floquet solution
        ss = obj.floquet(A_0, A_c, A_s, omega);
        obj.steady_state_internal = ss;
      else                                       % If the above odes are time independent, 
        assert(obj.is_stable, "There is no steady state covariance matrix, as the system is not stable!"); % Check if there exist a steady state!
        
        R_ss = linsolve(obj.A, -obj.N);          % Calculate steady-state mean quadratures
        V_ss = lyap(obj.A, obj.D);               % Calculate steady-state covariance matrix
        
        obj.steady_state_internal = gaussian_state(R_ss, V_ss); % Generate the steady state
        ss = obj.steady_state_internal;          %
      end
    end
    
    function ss = floquet(obj, A_0, A_c, A_s, omega)
      % Calculates the staeady state of a system with periodic Hamiltonin/drift matrix
      % Uses first order approximation in Floquet space for this calculation
      %
      % Higher order approximations will be implemented in the future
      % 
      % PARAMETERS:
      %   obj   - class instance
      % 
      %   A_0, A_c, A_s - components of the Floquet decomposition of the drift matrix
      %   omega - Frequency of the drift matrix
      % 
      % CALCULATES:
      %   obj.steady_state_internal with the steady state (gaussian_state)
      %   ss - gaussian_state with steady state of the system
      
      M = obj.Size_matrices;                     % Size of the time-dependent matrix
      Id = eye(M);                               % Identity matrix for the system size
      
      A_F = [[A_0,    A_c   ,     A_s  ];
             [A_c,    A_0   , -omega*Id];
             [A_s, +omega*Id,     A_0  ]];       % Floquet drift     matrix
      
      D_F = blkdiag(obj.D, obj.D, obj.D);        % Floquet diffusion matrix
      
      N_F = repmat(obj.N, 3, 1);                 % Floquet mean noise vector
      
      R_ss_F = linsolve(A_F, -N_F);              % Calculate steady-state Floquet mean quadratures vector
      V_ss_F = lyap(A_F, D_F);                   % Calculate steady-state Floquet covariance matrix
      
      R_ss = R_ss_F(1:M);                        % Get only the first entries
      V_ss = V_ss_F(1:M, 1:M);                   % Get only the first sub-matrix
      
      obj.steady_state_internal = gaussian_state(R_ss, V_ss); % Generate the steady state
      ss = obj.steady_state_internal; 
    end
    
    function langevin_semi_classical(obj, t_span, N_ensemble)
      % Solve the semi-classical Langevin equation for the expectation value of the quadrature operators
      % using a Monte Carlos simulation to numericaly integrate the Langevin equations
      % 
      % The initial conditions follows the initial state probability density in phase space
      % The differential stochastic equations are solved through a Euler-Maruyama method
      %
      % PARAMETERS:
      %   obj   - class instance
      %   N_ensemble (optional) - number of iterations for Monte Carlos simulation, default value: 200
      %
      % CALCULATES:
      %   obj.R_semi_classical - matrix with the quadratures expectation values of the time evolved system where 
      %   obj.R_semi_classical(i,j) is the i-th quadrature expectation value at the j-th time
      
      obj.t = t_span;                                      % Timestamps for the simulation
      obj.N_time = length(t_span);                         % Number of timestamps
      
      if nargin < 3
        N_ensemble = 2e+2;                                 % Number of Monte Carlos iterations
      end
      
      dt = obj.t(2) - obj.t(1);                            % Time step
      sq_dt = sqrt(dt);                                    % Square root of time step (for Wiener proccess in the stochastic integration)
      
      noise_amplitude = obj.N + sqrt(diag(obj.D));         % Amplitude for the noises (square root of the auto correlations)
      
      mean_0 = obj.initial_state.R;                        % Initial mean value
      std_deviation_0 = sqrt(diag(obj.initial_state.V));   % Initial standard deviation
      
      obj.R_semi_classical = zeros([obj.Size_matrices, obj.N_time]); % Matrix to store each quadrature ensemble average at each time
      
      if isa(obj.A,'function_handle')                      % I have to check if there is a time_dependency on the odes
        AA = @(t) obj.A(t);                                % Rename the function that calculates the drift matrix at each time
      else
        AA = @(t) obj.A;                                   % If A is does not vary in time, the new function always returns the same value 
      end
      
      for i=1:N_ensemble                                   % Loop on the random initial positions (% Monte Carlos simulation using Euler-Maruyama method in each iteration)
        
        X = zeros([obj.Size_matrices, obj.N_time]);        % For this iteration, this matrix stores each quadrature at each time (first and second dimensions, respectively)
        X(:,1) = normrnd(mean_0, std_deviation_0);         % Initial Cavity position quadrature (normal distribution for vacuum state)
        
        noise = randn(size(X));
        for k = 1:obj.N_time-1                             % Euler-Maruyama method for stochastic integration
          X(:,k+1) = X(:,k) + AA(obj.t(k))*X(:,k)*dt + sq_dt*noise_amplitude.*noise(:,k);
        end                                                % loop on time
        
        obj.R_semi_classical = obj.R_semi_classical + X;   % Add the new  Monte Carlos iteration quadratures to the same matrix
      end                                                  % loop on each Monte Carlos realization
      
      obj.R_semi_classical = obj.R_semi_classical/N_ensemble; % Divide the ensemble sum to obtain the average quadratures at each time
      
      % fprintf("Langevin simulation ended\n")             % Warn user that heavy calculation started
    end
    
    %%%%%%%%%%% Plotting functions, will not be present in final version %%%%%%%%%
    
    function plot_semi_classical(obj)
       mode_colors = [[0, 0, 0] ; [0, 0, 1]; [1, 0, 0]; [0, 0.5, 0]; [0.75, 0, 0.75]; [0.75, 0.75, 0]; [0.6350, 0.0780, 0.1840]]; % Array with colors to be consistently used throughout the plots
       
       if isempty(obj.R_semi_classical)
        disp("No semi-classical quadratures calulated, please calculate it before plotting it!")
        return
      end
      
      N_particles = obj.Size_matrices/2.0 - 1.0;
      
      figure('Name', "Semi-classical Quadratures Phase Space");           % Open a figure to contain the subplots
      
      x_label = "x_" + string(1:N_particles);          % Create the labels for the x axis
      x_label = ["q", x_label];
      
      y_label = "p_" + string(1:N_particles);          % Create the labels for the y axis
      y_label = ["p", y_label];
      
      title_name = "Particle " + string(1:N_particles);% Create the title (each mode name)
      title_name = ["Cavity", title_name];
      
      for i=0:N_particles                              % For each mode
        x = obj.R_semi_classical(2*i+1, :);            % Position quadrature for the i-th mode (semi-classical)
        p = obj.R_semi_classical(2*i+2, :);            % Momentum quadrature for the i-th mode (semi-classical)
        
        dq = x(2:end)-x(1:end-1);                          % Increments in the position for each time step to be used in the arrow plot
        dp = p(2:end)-p(1:end-1);                          % Increments in the momentum for each time step to be used in the arrow plot
        
        quiver_scale = 0;                                  % Tell MATLAB to not scale the arrows
        
        subplot(1, N_particles+1, i+1);                % Specify why subplot to use
        quiver(x(1:end-1), p(1:end-1), dq, dp, quiver_scale, 'Color', mode_colors(i+1,:)) % Arrow plot the phase space trajectory
        hold on
        plot(x(1), p(1), 'x', 'Color', mode_colors(i+1,:), 'MarkerSize', 12)
        
        title(title_name(i+1))                             % Add a title
        xlabel(x_label(i+1));                              % Label the x axis
        ylabel(y_label(i+1),'rotation',0,'VerticalAlignment','middle'); % Label the p axis
        
        set(gca, 'Fontsize', 18)                           % Set text to be on proper font size
        
      end
      
    end
    
    function plot_mean_quadratures(obj)
       mode_colors = [[0, 0, 0] ; [0, 0, 1]; [1, 0, 0]; [0, 0.5, 0]; [0.75, 0, 0.75]; [0.75, 0.75, 0]; [0.6350, 0.0780, 0.1840]]; % Array with colors to be consistently used throughout the plots
       
       if isempty(obj.R)
        disp("No mean quadratures calulated, please calculate it before plotting it!")
        return
      end
      
      N_particles = obj.Size_matrices/2.0 - 1.0;
      
      figure('Name', "Mean Quadratures Phase Space");           % Open a figure to contain the subplots
      
      x_label = "x_" + string(1:N_particles);          % Create the labels for the x axis
      x_label = ["q", x_label];
      
      y_label = "p_" + string(1:N_particles);          % Create the labels for the y axis
      y_label = ["p", y_label];
      
      title_name = "Particle " + string(1:N_particles);% Create the title (each mode name)
      title_name = ["Cavity", title_name];
      
      for i=0:N_particles                              % For each mode
        x = obj.R(2*i+1, :);                           % Position quadrature for the i-th mode
        p = obj.R(2*i+2, :);                           % Momentum quadrature for the i-th mode
        
%         x_semi = obj.R_semi_classical(2*i+1, :);       % Position quadrature for the i-th mode (semi-classical)
%         p_semi = obj.R_semi_classical(2*i+2, :);       % Momentum quadrature for the i-th mode (semi-classical)
        
        subplot(1, N_particles+1, i+1);                % Specify why subplot to use
        
        %hold on
        %plot(x(1), p(1), 'x', 'Color', mode_colors(i+1,:), 'MarkerSize', 12)
        
        if all(x == x(1)) && all(p == p(1)) % This happens for initial
          plot(x(1), p(1), 'x', 'Color', mode_colors(i+1,:), 'MarkerSize', 12)
        else
          plot(x, p, 'Color', mode_colors(i+1,:))        % Arrow plot the phase space trajectory
          hold on
          
%           plot(x, p, '--', 'Color', mode_colors(i+1,:))        % Arrow plot the phase space trajectory
%           hold on
%           plot(x_semi, p_semi, 'Color', mode_colors(i+1,:))
          
          plot(x(1), p(1), 'x', 'Color', mode_colors(i+1,:), 'MarkerSize', 12)
        end
        
        title(title_name(i+1))                             % Add a title
        xlabel(x_label(i+1));                              % Label the x axis
        ylabel(y_label(i+1),'rotation',0,'VerticalAlignment','middle'); % Label the p axis
        
        set(gca, 'Fontsize', 18)                           % Set text to be on proper font size
        
      end
      
    end
    
    function plot_CM(obj, skip)
      % Visualize the time evolution of the covariance matrix 
      % This is a very bad way of understanding the dynamics,
      % but can be used to see if the code is actually performing any time evolution
      
      for i=1:skip:size(obj.V)
        heatmap(obj.V{i});
        drawnow
      end
    end
    
  end % end methods
end % end class



% Auxiliar function defining the Lyapunov equation
function dVdt_vector = lyapunov_ode(~, V_old_vector, A, D)

M = size(A, 1);                            % System dimension (N_particles + 1 cavity field)partículas  + 1 campo)

V_old = reshape(V_old_vector, [M, M]);     % Vector -> matrix

dVdt = A*V_old + V_old*transpose(A) + D;   % Calculate how much the CM derivative in this time step
 
dVdt_vector = reshape(dVdt, [M^2, 1]);     % Matrix -> vector

end







% DRAFT:

%       for i=1:n
%         A_F = blkdiag(A_F, A_0);
%       end
%       ss = 1;
      

%       V_ss_numerical = obj.V{end};             % Find the last CM found numerically
%       V_diff = V_ss - V_ss_numerical;
%       relative_error = mean(abs(V_diff), 'all')./mean(abs(V_ss), 'all');
%       
%       if max(relative_error) < 1e-2            % Test if the system achieved the steady-state and warn user !
%         fprintf("Steady state achieved!\n")
%       else
%         fprintf("Steady state has not be achieved!\n")
%       end



% function nbar = occupation_number(obj)
% nbar = zeros(obj.Size_matrices/2, obj.N_time);
% 
% for i=1:length(obj.V)
%   Variances = diag(obj.V{i});              % From the current CM, take take the variances
%   
%   mean_x = obj.R(1:2:end, i);               % Odd  entries are the mean values of the position
%   mean_p = obj.R(2:2:end, i);               % Even entries are the mean values of the momentum
%   
%   Var_x = Variances(1:2:end);              % Odd  entries are position variances
%   Var_p = Variances(2:2:end);              % Even entries are momentum variances
%   
%   nbar(:, i) = 0.25*( Var_x + mean_x.^2 + Var_p + mean_p.^2 ) - 0.5; % Calculate occupantion numbers at current time
% end
% end


%         temp2(i).R = RR(:, i);
%         temp2(i).V = VV{i};
        
%         assert(all(temp2(i).R == temp(i).R), "R!")
%         assert(all(temp2(i).V == temp(i).V, 'all'), "V!")
%         assert(all(temp2(i).Omega == temp(i).Omega, 'all'), "Omega!")
%         assert(all(temp2(i).N_modes == temp(i).N_modes), "N_modes!")


% obj.t = t_span;                             % Store the information about the time interval for the simulation
%       obj.N_time = length(t_span);                % and its length on the simulation object


