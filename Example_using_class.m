%% Parameters adapted from Cavity Cooling of a Levitated Nanosphere by Coherent Scattering
omega =    2*pi*[190e+3]%;  160e+3;  180e+3];       % Natural frequency particle 1        [Hz]
g     =    2*pi*[ 42e+3]%; 35.3e+3; 19.8e+3];       % Coupling strength                   [Hz]
gamma =    2*pi*[ 10e+3]%;   10e+3;   10e+3];       % Damping                             [Hz]
T     =         [  300]%;    1e-6;    1e-6];       % Initial temperature of each particle [K]
T_environment = [  330]%;    1e-6;    1e-6];       % Temperature for the environment of each particle [K]

Delta = +omega(1);                       % Cavity-tweezer detuning             [Hz]
kappa = 2*pi*193e+2;                     % Cavity linewidth                    [Hz]

% Para só considerar a dinamica unitária, faça kappa = gamma = 0 !

%% Properties for the time interval under study
[T_idx, idx] = max(T);
N = 30;
t = linspace(0, N/omega(idx), 1e+4);     % Time stamps (for N oscilations)


%% Example of simulation
many_particles = simulation(omega, g, gamma, T, T_environment, Delta, kappa(end)); % Create a simulation variable
many_particles.run(t, "occupation_number", "steady_state");                 % Run the created simulation

% You can choose to calculate only what suits you through extra accepted parameters:
%["occupation_number", "heat_flux", "entanglement", "entropy", "steady_state", "langevin", "fidelity_test"]

many_particles.plot();                  % Plot the results for N particles

% You can also choose to plot things individually
% plot_entanglement , plot_entanglement_and_entropy , plot_single_mode_entropy , plot_occupation_number , plot_heat_fluxes , plot_phase_space, plot_fidelity_approx


%% Compare with N runs with a single particle
% hold on
% 
% for idx = 1:many_particles.N_particles
%   single_particle = simulation(omega(idx), g(idx), gamma(idx), T(idx), T_environment(idx), Delta, kappa);
%   single_particle.run(t, "occupation");
%   
%   plot(single_particle.t, single_particle.particles{1}.nbar, '--', 'Linewidth', 1.5, 'DisplayName', "Single Particle " + idx, 'Color', many_particles.mode_colors(idx+1,:))
% end



