%% Add path to auxiliar functions (these are important for CS simulation, but not for the toolbox)
addpath('auxiliar')


%% Universal Constants
k_B = 1.380649e-23;                      % Boltzmann's constant      [J/K]
hbar = 1.054571817e-34;                  % Reduced Planck's constant [J*s]


%% Experimental parameters (adapted from Aspelmeyer's ground state cooling article)
R = 100e-9;                                 % Particle radius [m]
P_t = 1000e-3;                              % Tweezer power   [W]
W_0 = 0.852e-6;                             % Cavity waist    [m]
T_gas = 130;                                % Environmental gas temperature [K]
nbar_0 = 0;%0.43/10;                        % Initial occupation number

% Particle natural frequency, coupling strength, damping, cavity-tweezer detunning, cavity linewidth, recoil heating rate and gas collision heating rate [Hz]
[omega_0, g_0, gamma_0, Delta_0, kappa_0, Gamma_recoil_0, Gamma_gas_0] = parameters(R, P_t, W_0, T_gas); % [Hz]

T_0 = hbar*omega_0/(k_B*log( (nbar_0+1)/nbar_0)); % Temperature for desired occupation number [K]

Gamma = Gamma_gas_0 + Gamma_recoil_0;       % Decoherence rate [Hz]
t_end = 1/Gamma;                            % Coherence time [s]


%% Build the matrices that define the dynamics of the system
N_particles = 2;

omega = omega_0*ones(N_particles, 1);       % Natural frequency of the particles    [Hz]
g     =     g_0*ones(N_particles, 1);       % Coupling strength                     [Hz]
gamma = gamma_0*ones(N_particles, 1);       % Damping                               [Hz]
T     =     T_0*ones(N_particles, 1);       % Initial temperature of each particle   [K]
T_environment = T_gas*ones(N_particles, 1); % Temperature for each environmental gas [K]

Delta = Delta_0;                            % Cavity-tweezer detuning               [Hz]
kappa = kappa_0;                            % Cavity linewidth                      [Hz]

 [A, D, N] = CS_simple_dynamics(omega, g, gamma, T_environment, Delta, kappa);
%[A, D, N] = CS_new_dynamics   (omega, g, gamma, T_environment, Delta, kappa);

%% Initial state using gaussian_state class
initial_state = gaussian_state("coherent", 1e6);
nbar = nbar_0*ones(N_particles, 1);       % Natural frequency of the particles    [Hz]
for i=1:N_particles
  thermal = gaussian_state("thermal", nbar(i));
  initial_state = initial_state.tensor_product(thermal);
end


%% Simulate the time evolution of the system
complete = time_evolution(A, D, N, initial_state); % Create a simulation variable for a system with 3 particles and an optical cavity

t = linspace(0, 3*t_end, 3e+3);                      % Time stamps for the simulation
complete.run(t);                                   % Run the created simulation (Optional parameters: ["langevin", "lyapunov", "semi_classical", "steady_state", "Slartibartfast"])

% complete.plot_mean_quadratures();

complete.langevin_semi_classical(t);
% complete.plot_semi_classical();

ss = complete.steady_state();





%% DRAFT

% Plot time evolving logarithmic negativity for mechanical bipartition
% log_neg = zeros(size(complete.t));
% for i=1:complete.N_time
% log_neg(i) = complete.state(i).logarithmic_negativity([2,3]);
% end
% plot(t, log_neg)

% Plot time evolution of covariance matrix 
% for i=1:size(complete.V)
% imshow(complete.V{i});
% drawnow
% end

% Decreasing fidelity
% base_state = gaussian_state("thermal", 0);
% f = zeros(1e3, 1);
% for i=0:1e3
%   target_state = gaussian_state("thermal", i);
%   f(i+1) = base_state.fidelity(target_state);
% end
% loglog(f)




%% OLD CODE
% Coherent state
% initial_state = gaussian_state("coherent", 2 + 1*1i);

% Occupation number
% nbar = 1./( exp(hbar*omega./(k_B*T)) - 1 );

% Time-dependent diffusion matrix
% AA = @(t) A*sin(omega(1)*t);

% for i=1:length(T)
%   thermal = gaussian_state("thermal", nbar(i));
%   initial_state = initial_state.tensor_product(thermal);
% end

% %% Input parameters
% omega =    2*pi*[305.4e+3;  305.4e+3;  305.4e+3];     % Natural frequency of the particles [Hz]
% g     =    2*pi*[ 64.0e+3;   93.2e+3;  109.2e+3];     % Coupling strength                  [Hz]
% gamma =    2*pi*[ 9.57e-4;   9.57e-4;   9.57e-4];     % Damping                            [Hz]
% T     =         [ 4.6e-6;    4.6e-6;     4.6e-6];     % Initial temperature of each particle             [K]
% T_environment = [   300 ;       300;        300];     % Temperature for the environment of each particle [K]
% 
% % The length of the vectors with the parameters for the particles define the number of particles in the simulation
% 
% Delta = 2*pi*315e+3;                                  % Cavity field natural frequency      [Hz]   (Cavity-tweezer detuning)
% kappa = 2*pi*193e+3;                                  % Cavity linewidth                    [Hz]

% If you wish to consider only the closed unitary dynamics of the system, just make gamma = 0 and kappa = 0 by uncommenting the following lines
% kappa=0
% gamma = zeros(size(omega));

% nbar = 1./( exp(hbar*omega./(k_B*T)) - 1 );

% %% Time interval for simulation
% t = linspace(0, 4.2e-6, 1e+3);                        % Time interval for the simulation   [s]
