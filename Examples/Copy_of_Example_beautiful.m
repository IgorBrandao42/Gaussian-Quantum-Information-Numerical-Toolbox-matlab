%% Universal Constants
k_B = 1.380649e-23;                              % Boltzmann's constant      [J/K]
hbar = 1.054571817e-34;                          % Reduced Planck's constant [J*s]

%% Parameters
% Aspelmeyer
% omega = 2*pi*305e+3;                           % Particle natural frequency [Hz]
% gamma = 2*pi*6.2998e-4;  %1e-6 mbar pressure   % Damping constant [Hz]

% Quidant
omega = 2*pi*197e+3;                             % Particle natural frequency [Hz]
gamma = 2*pi*881.9730;     % 1.4 mbar pressure   % Damping constant [Hz]

T_env = 300;                                     % Environmental temperature
T_0   = 300;                                       % Initial particle temperature

nbar_env = 1/( exp(hbar*omega/(k_B*T_env)) - 1 );% Environmental    occupation number
nbar_0   = 1/( exp(hbar*omega/(k_B*T_0  )) - 1 );% Initial particle occupation number



A =[[    0   ,  +omega ];                        % Drift matrix for harmonic potential
    [ -omega ,  -gamma ]];
        
D = diag([0, 2*gamma*(2*nbar_env+1)]);           % Diffusion matrix
N = zeros(2,1);                                  % Mean noise vector



%% Initial state
initial = gaussian_state("thermal", nbar_0);     % Initial state
% r = 1.1;
% initial.displace(300*1i-450);
% initial.squeeze(r);
% initial.rotate(-pi/4);
% nbar_0 = initial.occupation_number;

%% Timestamps for the dynamics
t = linspace(0, 2*pi/omega, 1e4);                % Timestamps for simulation


%% Simulation
simulation = time_evolution(A, D, N, initial);   % Simulate!
states = simulation.run(t);


%% Plot 
variance_W = 20*min([initial.V(1,1), initial.V(2,2)]);

x = sqrt(variance_W)*linspace(-1, +1, 150);                 % Region to plot wigner function
p = sqrt(variance_W)*linspace(-1, +1, 150);
[X, P] = meshgrid(x, p);

figure(2)
skip = 50;                                       % Make animation faster by increasing this!
for i=1:length(states)                      % Play animation  
  
  W = states(i).wigner(X, P);                    % Calculate the wigner function
  
  surf(X,P,W)                                    % Make pretty plot
  shading interp
  view(0,90)
  xlabel('x')
  xlim([x(1), x(end)])
  ylim([p(1), p(end)])
  ylabel('p', 'rotation', 0, 'VerticalAlignment', 'middle');
  axis square
  colorbar
  drawnow                                        % Display it
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
