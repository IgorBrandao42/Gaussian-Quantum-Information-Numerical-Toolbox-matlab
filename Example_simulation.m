%% Input parameters
omega =    2*pi*[305.4e+3;  305.4e+3;  305.4e+3];     % Natural frequency of the particles [Hz]
g     =    2*pi*[ 64.0e+3;   93.2e+3;  109.2e+3];     % Coupling strength                  [Hz]
gamma =    2*pi*[ 9.57e-4;   9.57e-4;   9.57e-4];     % Damping                            [Hz]
T     =         [ 4.6e-6;    4.6e-6;     4.6e-6];     % Initial temperature of each particle             [K]
T_environment = [   300 ;       300;        300];     % Temperature for the environment of each particle [K]

% The length of the vectors with the parameters for the particles define the number of particles in the simulation

Delta = 2*pi*315e+3;                                  % Cavity field natural frequency      [Hz]   (Cavity-tweezer detuning)
kappa = 2*pi*193e+3;                                  % Cavity linewidth                    [Hz]

% If you wish to consider only the closed unitary dynamics of the system, just make gamma = 0 and kappa = 0 by uncommenting the following lines
% kappa=0
% gamma = zeros(size(omega));


%% Time interval for simulation
t = linspace(0, 4.2e-6, 1e+3);                        % Time interval for the simulation   [s]


%% Example of a complete simulation
complete = simulation(omega, g, gamma, T, T_environment, Delta, kappa); % Create a simulation variable for a system with 3 particles and an optical cavity

complete.run(t);                                 % Run the created simulation

complete.plot();                                 % Plot the results


%% Example of a simple simulation with different number of particles (uncomment to try it out!)
% simple = simulation(omega(1:2), g(1:2), gamma(1:2), T(1:2), T_environment(1:2), Delta, kappa); % Simulation with the first two particles
% 
% simple.run(t, "occupation_number", "langevin");% Run only the simulation of occupation number and quadratures
% 
% simple.plot_occupation_number();               % Plot the time evolving occupation number for each mode
% simple.plot_phase_space();                     % Plot the time evolving expectation values of the quadrature on the phase space


%% Information about the cavity and particles can be easily recovered (uncomment to try it out!)
complete.cavity                                  % Optical cavity of the complete simulation
complete.particles{1}                            % Particle 1     of the complete simulation
% 
% simple.cavity                                  % Optical cavity of the simple simulation
% simple.particles{1}                            % Particle 1     of the simple simulation



