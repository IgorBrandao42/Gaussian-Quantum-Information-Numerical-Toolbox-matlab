function [W_particles, W_cavity] = phonon_work_fluxes_RWA(states, g, delta, t)
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
    W_particles(1, i) = - 0.5*g(1)*( states(i).V(1, 2*1+2) + states(i).R(1)*states(i).R(2*1+2) );
    
    W_particles(2, i) = - 0.5*g(2)*( states(i).V(1, 2*2+2) + states(i).R(1)*states(i).R(2*2+2) )*cos( delta*t(i) )...
                        + 0.5*g(2)*( states(i).V(2, 2*2+2) + states(i).R(2)*states(i).R(2*2+2) )*sin( delta*t(i) );
    
    W_cavity(i)       = - 0.5*g(1)*( states(i).V(2, 2*1+1) + states(i).R(2)*states(i).R(2*1+1) )...
                        - 0.5*g(2)*( states(i).V(2, 2*2+1) + states(i).R(2)*states(i).R(2*2+1) )*cos( delta*t(i) )...
                        - 0.5*g(2)*( states(i).V(1, 2*2+1) + states(i).R(1)*states(i).R(2*2+1) )*sin( delta*t(i) );
end

end



