function [W_particles, W_cavity] = phonon_work_fluxes(states, g, nu, t)
% Calculates the heat fluxes from the environments into each mode
%
% INPUTS
%    states   - time evolved gaussian states
%    kappa    - cavity linewidth
%    gamma    - damping coefficient
%    nbar_env - occupation number for the environment of each particle
%
% RETURNS:
%    W_particles(j, i) - phonon heat flux from the environment to the j-th particle at the i-th timestamp
%    W_cavity(i)       - phonon heat flux from the environment to the cavity mode

N_particles     = states(1).N_modes - 1;                   % The number of particles is the total number of modes minus the number of cavity modes (1)

W_cavity    = zeros(size(t));                              % Initialize variable to store heat flux from environment to the  cavity
W_particles = zeros(N_particles, length(t));               % Initialize variable to store heat flux from environment to each particle 

for i=1:length(states)                                     % Loop trough each time evolved state and calculate eac heat flux
  for j = 1:N_particles                                    % Loop through each particle
    W_particles(j, i) =       -2.0*g(j)*( states(i).V(1, 2*j+2) + states(i).R(1)*states(i).R(2*j+2) )*cos( nu(j)*t(i) );
    
    W_cavity(i) = W_cavity(i) -2.0*g(j)*( states(i).V(2, 2*j+1) + states(i).R(2)*states(i).R(2*j+1) )*cos( nu(j)*t(i) );
  end
end

end



