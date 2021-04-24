function [A, D, N] = CS_time_dependency(omega, g, gamma, T_environment, Delta, kappa)
% Calculates the parameters for the Langevin and Lyapunov equations dictating the dynamics
% of an experiment with N nanoparticles interacting with an optical cavity field
% See https://arxiv.org/abs/2102.08969
%
% PARAMETERS:
%   omega         - Array with Natural frequencies  of each nanoparticle
%   g             - Array with Coupling strength   for each nanoparticle's interaction with the cavity field
%   gamma         - Array with Mechanical damping  for each nanoparticle
%   T_init        - Array with Initial temperature for each particle
%   T_environment - Array with Temperature for each heat bath
%   Delta         - Natural frequency of the optical mode (in the context of Coherent Scattering, this is the cavity-tweezer detunning)
%   kappa         - Cavity linewidth
%
% The number of particles is inferred from the size of the input parameters (they must all be the same!)
%
% CALCULATES:
%   A      - the diffusion matrix of the system
%   D      - the drift matrix of the system
%   N      - the initial covariance matrix
%   stable - boolean telling if the system is stable or not

%% Make sure the user gave the right number of parameters
number_parameters = [length(omega), length(g), length(gamma), length(T_environment)];
assert(all(number_parameters == number_parameters(1)), "Wrong size of vectors with parameters! Impossible to determine number of particles, ending program...")

%% Universal Constants
k_B = 1.380649e-23;                              % Boltzmann's constant      [J/K]
hbar = 1.054571817e-34;                          % Reduced Planck's constant [J*s]

%% Useful constants 
N_particles = length(omega);                     % Number of particles
Size_matrices = 2*(N_particles+1);               % Size of system

nbar_environment = 1./( exp(hbar*omega./(k_B*T_environment)) - 1 ); % Occupation number for the environment (heat bath)


%% Drift matrix
A = [[-kappa/2 , Delta]; [-Delta   , -kappa/2]]; % Matrix defining the Langevin equation (without noises)
for i=1:N_particles
  A_particle = [[0, omega(i)]; [-omega(i), -gamma(i)]];
  A = blkdiag(A, A_particle);
  A(2      , 2*i+1) = -2*g(i);
  A(2*(i+1),   1    ) = -2*g(i);
end


%% Diffusion matrix
D_diag = zeros([Size_matrices, 1]);              % Diagonal matrix with amplitude of noises' autocorrelators
D_diag(1:2) = [kappa; kappa];
D_diag(4:2:end) = 2*gamma.*(2*nbar_environment + 1);
D = diag(D_diag);


%% Mean value of the noises
N = zeros(Size_matrices, 1);                     % Mean of the noises acting on the system (will probabliy always be null!)

end