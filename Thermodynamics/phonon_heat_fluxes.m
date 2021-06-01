function [J_env_particles, J_env_cavity] = phonon_heat_fluxes(states, kappa, gamma, nbar_env)
% Calculates the heat fluxes from the environments into each mode
%
% INPUTS
%    states   - time evolved gaussian states
%    kappa    - cavity linewidth
%    gamma    - damping coefficient
%    nbar_env - occupation number for the environment of each particle
%
% RETURNS:
%    J_env_particles(j, i) - phonon heat flux from the environment to the j-th particle at the i-th timestamp
%    J_env_cavity(i)       - phonon heat flux from the environment to the cavity mode

N_particles     = states(1).N_modes - 1;                   % The number of particles is the total number of modes minus the number of cavity modes (1)

J_env_cavity    = zeros(size(states));                     % Initialize variable to store heat flux from environment to the  cavity
J_env_particles = zeros(N_particles, length(states));      % Initialize variable to store heat flux from environment to each particle 

for i=1:length(states)                                     % Loop trough each time evolved state and calculate eac heat flux
  J_env_cavity(i) = kappa/4*( 1 - states(i).V(1, 1) ) + kappa/4*( 1 - states(i).V(2, 2) );   
  
  for j = 1:N_particles                                    % Loop through each particle
    J_env_particles(j, i) = gamma(j)/2*( 2*nbar_env(j) + 1 - states(i).V(2*j+2, 2*j+2) );
  end
end

end









% N_particles = states(1).N_modes - 1;                       % The number of particles is the total number of modes minus the number of cavity modes (1)
% J_env_cav_Q = zeros(size(t));                              % Initialize variable to store heat flux from environment to position quadrature of the cavity
% J_env_cav_P = zeros(size(t));                              % Initialize variable to store heat flux from environment to momentum quadrature of the cavity
% J_env       = zeros(N_particles, length(t));               % Initialize variable to store heat flux from environment to each particle 
% 
% for i=1:length(states)                                     % Loop trough each time evolved state and calculate eac heat flux
%   J_env_cav_Q(i) = kappa/4*( 1 - states(i).V(1, 1) );   
%   J_env_cav_P(i) = kappa/4*( 1 - states(i).V(2, 2) );
%   
%   for j = 1:N_particles
%     J_env(j, i) = gamma(j)/2*( 2*nbar_env(j) + 1 - states(i).V(2*j+2, 2*j+2) );
%   end
% end
