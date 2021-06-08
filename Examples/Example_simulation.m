%% Universal Constants
k_B  = 1.380649e-23;                        % Boltzmann's constant      [J/K]
hbar = 1.054571817e-34;                     % Reduced Planck's constant [J*s]
c    = 299792458;                           % Speed of light [m/s]

%% Experimental parameters (adapted from Aspelmeyer's ground state cooling article)
R = 100e-9;                                 % Particle radius [m]
P_t = 1000e-3;                              % Tweezer power   [W]
W_0 = 0.852e-6;                             % Cavity waist    [m]
T_gas = 300;                                % Environmental gas temperature [K]

% Particle natural frequency, coupling strength, damping, cavity-tweezer detunning, cavity linewidth, recoil heating rate and gas collision heating rate [Hz]
[omega_0, g_0, gamma_0, Delta_0, kappa_0, Gamma_recoil_0, Gamma_gas_0, E_d] = parameters(R, P_t, W_0, T_gas); % [Hz]

% gamma_0 = 0;
% kappa_0 = 5*kappa_0;

% T_0 = hbar*omega_0/(k_B*log( (nbar_0+1)/nbar_0)); % Temperature for desired occupation number [K]

Gamma = Gamma_gas_0 + Gamma_recoil_0;       % Decoherence rate [Hz]
t_end = 25/Gamma;                            % Coherence time [s]


%% Build the matrices that define the dynamics of the system
N_particles = 2;

omega = omega_0*ones(N_particles, 1);       % Natural frequency of the particles    [Hz]
g     =     g_0*ones(N_particles, 1);       % Coupling strength                     [Hz]
gamma = gamma_0*ones(N_particles, 1);       % Damping                               [Hz]
% T     =                [T_0; 10*T_0];       % Initial temperature of each particle   [K]
T_environment = T_gas*ones(N_particles, 1); % Temperature for each environmental gas [K]

Delta = Delta_0;                            % Cavity-tweezer detuning               [Hz]
kappa = kappa_0;                            % Cavity linewidth                      [Hz]

delta = sum(omega);                         % Tweezer-tweezer detuning [Hz]
omega_cav = c/780e-9*2*pi;                  % Cavity frequency [Hz]
nu = [omega_cav - omega(1); omega_cav + omega(2)]; % Tweezer frequencies [Hz]

% [A, D, N] = CS_simple_dynamics(omega, g, gamma, T_environment, Delta, kappa);
% [A, D, N] = CS_time_dependency(omega, g, gamma, T_environment, omega_cav, kappa, nu);
[A, D, N] = CS_time_dependency_RWA(omega, g, gamma, T_environment, Delta, kappa, delta);


%% Initial state using gaussian_state class
initial_state = gaussian_state("coherent", 0);
T_0 = T_gas; % Particles' initial temperature [K]
nbar_0 = 1./( exp(hbar*omega/(k_B*T_0)) - 1 );                      % Initial occupation number           
for i=1:N_particles
  thermal = gaussian_state("thermal", nbar_0(i));
  initial_state = initial_state.tensor_product(thermal);
end

%% Simulate the time evolution of the system
complete = time_evolution(A, D, N, initial_state); % Create a simulation variable for a system with 3 particles and an optical cavity

              % 1000*2*pi/nu(1)
t = linspace(0, t_end, 3e+4);                       % Time stamps for the simulatio
complete.run(t);

states = complete.state;



%% Calculate some properties
% log_neg_mec_opto = zeros(complete.N_time,1);
% log_neg_opto_mec = zeros(complete.N_time,1);
% log_neg_mec = zeros(complete.N_time,1);
% % duan = zeros(complete.N_time,1);
% 
% for i=1:complete.N_time
%   log_neg_opto_mec(i) = states(i).logarithmic_negativity([1,2]);
%   log_neg_mec_opto(i) = states(i).logarithmic_negativity([1,3]);
%   log_neg_mec(i)      = states(i).logarithmic_negativity([2,3]);
% % duan(i)    = states(i).duan_criteria([2,3]);
% end

nbar_env = 1./( exp(hbar*omega/(k_B*T_gas)) - 1 );         % Environmental occupation number

nbar = zeros(complete.Size_matrices/2, complete.N_time);
for i=1:complete.N_time
nbar(:,i) = states(i).occupation_number();
end


%% Plot those properties 
mode_colors = [[0, 0, 0] ; [0, 0, 1]; [1, 0, 0]; [0, 0.5, 0]; [0.75, 0, 0.75]; [0.75, 0.75, 0]; [0.6350, 0.0780, 0.1840]]; % Array with colors to be consistently used throughout the plots
mode_name = "Particle " + string(1:N_particles); % Create the begining of the title for each mode
mode_name = ["Cavity", mode_name];

figure(2)
clf 
hold on
for i=1:complete.Size_matrices/2
  plot(t, nbar(i, :), 'Linewidth', 1.5, 'DisplayName', mode_name(i), 'Color', mode_colors(i,:))
end
ylabel('$\bar{n}$','Interpreter','Latex','rotation',0,'VerticalAlignment','middle');
xlabel('t [s]');
xlim([t(1), t(end)]);

legend
set(gca, 'Fontsize', 14)
set(gca, 'TickDir', 'out')
title("Occupation number")
% set(gca, 'YScale', 'log')






%% Calculate and plot phonon heat flux from the environments

heat_name = ["J_{env, 1}", "J_{env, 2}"];

[J_env_particles, J_env_cavity] = phonon_heat_fluxes(states, kappa, gamma, nbar_env);

figure(3); clf
plot(t, J_env_cavity, 'DisplayName', "J_{env, cav}", 'Color', mode_colors(1,:), 'Linewidth', 1)
title("Phonon heat flux from the environment to the cavity")
legend

figure(4); clf; hold on
for j = 1:N_particles
  plot(t, J_env_particles(j, :), 'DisplayName', heat_name(j), 'Color', mode_colors(j+1,:), 'Linewidth', 1)
end
yline(0, 'k--', 'Linewidth', 2);
title("Phonon heat flux from the environment to each particle")
legend

%% Calculate and plot phonon heat flux from the environments (in the non-rotating frame)
work_name = ["W_{1}", "W_{2}"];

[W_particles, W_cavity] = phonon_work_fluxes_RWA(states, g, delta, t);
%[W_particles, W_cavity] = phonon_work_fluxes(states, g, nu, t);

figure(5); clf
plot(t, W_cavity, 'DisplayName', "W_{cav}", 'Color', mode_colors(1,:), 'Linewidth', 1)
title("Phonon work flux from particles to cavity  (non-rotating frame)")
legend

figure(6); clf; hold on
for j = 1:N_particles
  plot(t, W_particles(j, :), 'DisplayName', work_name(j), 'Color', mode_colors(j+1,:), 'Linewidth', 1)
end
yline(0, 'k--', 'Linewidth', 2);
title("Phonon work flux from cavity to each particle  (non-rotating frame)")
legend





%% DRAFT - old code
% J_env_cav = zeros(2, length(t));
% for i=1:length(states)
%     J_env_cav(1, i) = kappa/4*( 1 - states(i).V(1, 1) );
%     J_env_cav(2, i) = kappa/4*( 1 - states(i).V(2, 2) );
% end
% plot(t, J_env_cav(1, :), 'DisplayName', "J_{env, cav}", 'Color', mode_colors(1,:) )
% plot(t, J_env_cav(2, :), 'DisplayName', "J_{env, cav}", 'Color', mode_colors(1,:) )


% Plot wigner function for single particle
% variance_W = 5*min([states(1).V(1,1), states(1).V(2,2)]);
% x = states(1).R(1) + 5*sqrt(variance_W)*linspace(-1, +1, 150);                 % Region to plot wigner function
% p = states(1).R(2) + 5*sqrt(variance_W)*linspace(-1, +1, 150);
% [X, P] = meshgrid(x, p);
% skip = 1;
% for i=1:skip:length(states)
% W = states(i).wigner(X,P);
% surf(X,P,W)                                    % Make pretty plot
% shading interp
% view(0,90)
% xlabel('x')
% xlim([x(1), x(end)])
% ylim([p(1), p(end)])
% ylabel('p', 'rotation', 0, 'VerticalAlignment', 'middle');
% axis square
% drawnow
% end



% is_max = islocalmax(log_neg_mec);
% t_max = t(is_max);
% 
% subplot(3,1,1)
% plot(t, log_neg_opto_mec, 'b', 'Linewidth', 2)
% for i=1:length(t_max)
%   xline(t_max(i), 'k--');
% end
% 
% subplot(3,1,2)
% plot(t, log_neg_mec_opto, 'b', 'Linewidth', 2)
% for i=1:length(t_max)
%   xline(t_max(i), 'k--');
% end
% 
% subplot(3,1,3)
% plot(t, log_neg_mec, 'b', 'Linewidth', 2)
% for i=1:length(t_max)
%   xline(t_max(i), 'k--'); 
% end


% complete.plot_mean_quadratures();

%complete.langevin_semi_classical(t);
% complete.plot_semi_classical();

%ss = complete.steady_state();





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
% nbar = 1./( exp(hbar*omega./(k_B*T)) - 1 ); % Occupation number

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