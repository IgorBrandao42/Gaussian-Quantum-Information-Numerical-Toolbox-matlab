%% Input parameters
omega =    2*pi*[190e+3;  160e+3;  180e+3];      % Natural frequency of the particles  [Hz]
g     =    2*pi*[ 42e+3; 35.3e+3; 39.8e+3];      % Coupling strength                   [Hz]
gamma =    2*pi*[ 10e+3;   10e+3;   10e+3];      % Damping                             [Hz]
T     =         [  1e-3;    1e-3;    1e-3];      % Initial temperature of each particle             [K]
T_environment = [  1e-3;    1e-3;    1e-3];      % Temperature for the environment of each particle [K]

Delta = +omega(1);                               % Cavity field natural frequency      [Hz]   (Cavity-tweezer detuning)
kappa = 2*pi*193e+3;                             % Cavity linewidth                    [Hz]

% The length of the vectors with the parameters for the particles define the number of particles in the simulation
% kappa=0;
% gamma = zeros(size(omega));
% If you wish to consider only the closed unitary dynamics of the system, just make gamma = 0 and kappa = 0


%% Time interval for simulation
t = linspace(0, 50/omega(1), 1e+5);                     % Time stamps for the calculations


%% Example of a complete simulation
% Create a simulation variable for a system with 3 particles and an optical cavity
complete = simulation(omega(1:2), g(1:2), gamma(1:2), T(1:2), T_environment(1:2), Delta, kappa);

complete.run(t, "entanglement", "mutual_information");  % Run the created simulation

%complete.plot_entanglement();
%complete.plot();                                 % Plot the results for N particles

% You can choose to calculate only what suits you through extra accepted parameters:
%["occupation_number", "lyapunov", "heat_flux", "entanglement", "mutual_information", "entropy", "steady_state", "langevin", "fidelity_test"]

% You can also choose to plot things individually
% plot_occupation_number, plot_mutual_information, plot_entanglement_and_entropy , plot_single_mode_entropy, plot_system_entropy , plot_entanglement , plot_heat_fluxes , plot_phase_space, plot_fidelity_approx


%% Example of a simple simulation with different number of particles
% simple = simulation(omega(1:2), g(1:2), gamma(1:2), T(1:2), T_environment(1:2), Delta, kappa); % Simulation with the first two particles
% 
% simple.run(t, "occupation_number", "langevin");  % Run only the simulation of occupation number and quadratures
% 
% simple.plot_occupation_number();                 % Plot the time evolving occupation number for each mode
% simple.plot_phase_space();                       % Plot the time evolving expectation values of the quadrature on the phase space
% % 
% % 
% % %% Information about the cavity and particles can be easily recovered
% complete.cavity                                % Optical cavity of the complete simulation
% complete.particles{1}                          % Particle 1     of the complete simulation
% 
% simple.cavity                                  % Optical cavity of the simple simulation
% simple.particles{1}                            % Particle 1     of the simple simulation



