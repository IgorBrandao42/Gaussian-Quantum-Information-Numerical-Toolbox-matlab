% Calculation of entanglement dynamics between two optically levitated nanspheres
% interacting through the Coherent-Scattering optomechanical interaction
% See https://arxiv.org/abs/2102.08969


%% Universal Constants
k_B  = 1.380649e-23;                        % Boltzmann's constant      [J/K]
hbar = 1.054571817e-34;                     % Reduced Planck's constant [J*s]


%% Parameters for the dynamics
omega = 2*pi*[420 ,  420]*1e3;              % Natural frequency of the particles    [Hz]
g     = 2*pi*[130 ,  130]*1e3;              % Coupling strength                     [Hz]
gamma = 2*pi*[9e-6, 9e-6];                  % Damping                               [Hz]
T_environment = [300, 300];                 % Temperature for each environmental gas [K]

Delta = 2*pi*315*1e3;                       % Cavity-tweezer detuning               [Hz]
kappa = 2*pi*193*1e3;                       % Cavity linewidth                      [Hz]


%% Build the matrices that define the dynamics of the system
% Drift matrix
A = [[-kappa/2 ,  Delta  ,     0     ,     0     ,     0     ,     0     ]; 
     [  -Delta , -kappa/2,  -2*g(1)  ,     0     ,  -2*g(2)  ,     0     ];
     [    0    ,    0    ,     0     , +omega(1) ,     0     ,     0     ];
     [ -2*g(1) ,    0    , -omega(1) , -gamma(1) ,     0     ,     0     ];
     [    0    ,    0    ,     0     ,     0     ,     0     , +omega(2) ];
     [ -2*g(2) ,    0    ,     0     ,     0     , -omega(2) , -gamma(2) ]];

% Diffusion matrix
nbar_env = 1./( exp(hbar*omega./(k_B*T_environment)) - 1 ); % Occupation number for the environment    
D_diag = zeros([6, 1]);                     % Diagonal of the diffusion matrix with the amplitude of noises' autocorrelators
D_diag(1:2) = [kappa; kappa];
D_diag(4:2:end) = 2*gamma.*(2*nbar_env + 1);
D = diag(D_diag);

% Mean value of the noises
N = zeros(6, 1);                                    % Mean of the noises acting on the system (will probabliy always be null!)


%% Initial state, simulated with the gaussian_state class
T_0 = 0;                                                 % Particles' initial temperature [K]
nbar_0 = 1./( exp(hbar*omega./(k_B*T_0)) - 1 );          % Initial occupation number    

initial_state = gaussian_state("vacuum");                % Initial state for the optical cavity
for i=1:2
  thermal = gaussian_state("thermal", nbar_0(i));        % Initial state for each particle
  initial_state = initial_state.tensor_product(thermal); % Initial global state is the tensor product of each single-mde state (no entanglement)
end


%% Simulate the time evolution of the system
simulation = gaussian_dynamics(A, D, N, initial_state); % Create a simulation variable for a system with 3 particles and an optical cavity

t = linspace(0, 1e-5, 1e+3);                % Time stamps for the simulatio
time_evolved_states = simulation.run(t);    % Calculate the time evolution as retrieve the time evolve gaussian states (array of gaussian_state instances)


%% Calculate mutual information, squeezing and entanglement
log_neg_mec_opto = zeros(simulation.N_time, 1);% Logarithmic negativity between cavity field and particle 1 for each timestamps
log_neg_opto_mec = zeros(simulation.N_time, 1);% Logarithmic negativity between cavity field and particle 2 for each timestamps
log_neg_mec      = zeros(simulation.N_time, 1);% Logarithmic negativity between the particles for each timestamps

I_particles = zeros(simulation.N_time, 1);     % Mutual information between the particles for each timestamps
I_total     = zeros(simulation.N_time, 1);     % Mutual information for the total state   for each timestamps

eta = zeros(simulation.N_time, 1);             % Squeezing parameter (ratio of squeezed ans anti-squeezed quadratures) for each timestamps

for i=1:simulation.N_time
  log_neg_opto_mec(i) = time_evolved_states(i).logarithmic_negativity([1,2]);
  log_neg_mec_opto(i) = time_evolved_states(i).logarithmic_negativity([1,3]);
  log_neg_mec(i)      = time_evolved_states(i).logarithmic_negativity([2,3]);
  
  I_total(i)      = time_evolved_states(i).mutual_information();
  
  mec_bipartition = time_evolved_states(i).partial_trace(1);
  I_particles(i)  = mec_bipartition.mutual_information();
  
  particle_1      =  time_evolved_states(i).only_modes(2);
  eta(i) = particle_1.squeezing_degree();
end


%% Plottting (there is no more usage of the toolbox below, only its results!)
% Colors and names for the following plot
mode_colors = [[0, 0, 0] ; [0, 0, 1]; [1, 0, 0]];   % Array with colors to be consistently used throughout the plots
mode_name = ["Cavity", "Particle 1", "Particle 2"]; % Name of each mode

% Plot logarithmic negativity
figure(1)
clf

subplot(3,1,1)
plot(t, log_neg_opto_mec, 'Linewidth', 1.5, 'DisplayName', "Log. Neg._{1,2}", 'Color', mode_colors(1,:))
legend
xline(t(end), 'k--', 'Linewidth', 1.5, 'DisplayName', 'Coherence time');
subplot(3,1,2)
plot(t, log_neg_mec_opto, 'Linewidth', 1.5, 'DisplayName', "Log. Neg._{1,3}", 'Color', mode_colors(2,:))
legend
xline(t(end), 'k--', 'Linewidth', 1.5, 'DisplayName', 'Coherence time');
subplot(3,1,3)
plot(t, log_neg_mec, 'Linewidth', 1.5, 'DisplayName', "Log. Neg._{2,3}", 'Color', mode_colors(3,:))
xline(t(end), 'k--', 'Linewidth', 1.5, 'DisplayName', 'Coherence time');
legend


% Plot squeezing parameter
figure(2)
clf
plot(t, eta, 'Linewidth', 1.5, 'DisplayName', "\eta(t)", 'Color', "k")
hold on
xline(t(end), 'k--', 'Linewidth', 1.5, 'DisplayName', 'Coherence time');
legend


figure(3)
subplot(2,1,1)
plot(t, I_total, 'Linewidth', 1.5, 'DisplayName', "I_{total}", 'Color', mode_colors(1,:))
xline(t(end), 'k--', 'Linewidth', 1.5, 'DisplayName', 'Coherence time');
legend

subplot(2,1,2)
plot(t, I_particles, 'Linewidth', 1.5, 'DisplayName', "I_{particles}", 'Color', mode_colors(3,:))
xline(t(end), 'k--', 'Linewidth', 1.5, 'DisplayName', 'Coherence time');
legend

