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
T_0   = 1;                                       % Initial particle temperature

nbar_env = 1/( exp(hbar*omega/(k_B*T_env)) - 1 );% Environmental    occupation number
nbar_0   = 1/( exp(hbar*omega/(k_B*T_0  )) - 1 );% Initial particle occupation number

%% Define the dynamics         
% A = @(t) [[                   0    ,  +omega*sin(8*omega*t) ]; % Drift matrix for harmonic potential
%           [ -omega*sin(8*omega*t) ,  -gamma ]];

A =[[    0   ,  +omega ];                        % Drift matrix for harmonic potential
    [ -omega ,  -gamma ]];
        
D = diag([0, 2*gamma*(2*nbar_env+1)]);           % Diffusion matrix
N = zeros(2,1);                                  % Mean noise vector

%% Initial state
initial = gaussian_state("thermal", nbar_0);     % Initial state
r = 3;
initial.squeeze(r);
initial.rotate(-pi/4);
nbar_0 = initial.occupation_number;

%% Timestamps for the dynamics
t = linspace(0, 2*pi/omega, 1e4);                   % Timestamps for simulation

%% Rough estimate of the total phonon lost
dt = t(end) - t(1);
J_env = gamma/2*(2*nbar_env+1 - initial.V(2,2));
phonon_lost_rough_estimate = J_env*dt

%% Simulation
simulation = time_evolution(A, D, N, initial);   % Simulate!
simulation.run(t);

states = simulation.state;                       % Retrieved time evolved states


%% Plot 
variance_W = 0.0003*min([initial.V(1,1), initial.V(2,2)]);

x = initial.R(1) + 200*sqrt(variance_W)*linspace(-1, +1, 250);                 % Region to plot wigner function
p = initial.R(2) + 200*sqrt(variance_W)*linspace(-1, +1, 250);
[X, P] = meshgrid(x, p);

figure(1)                                        % Figure for plotting
clf
nbar = zeros(size(t));                           % Variable to store occupation numbers

for i=1:1:length(states)
nbar(i) = states(i).occupation_number();
end
% 
h = zeros(1,3);
h(1) = semilogy(t, nbar, 'k', 'Linewidth', 2);
hold on
h(2) = yline(nbar_env, 'r--', 'linewidth', 3, 'DisplayName', 'Environment  occupation');
if nbar_0 ~= 0
  h(3) = yline(nbar_0  , 'b--', 'linewidth', 3, 'DisplayName', 'Particle initial occupation');
end
title('Time dependent occupation number')
xlabel('t [s]')
ax = gca;
ax.FontSize = 18;
legend(h(2:3))

% figure(4)
% clf
% 
% W2 = initial.wigner(X,P);
% surf(X, P, W2)
% view(0,90)
% shading interp
% xlim([x(1), x(end)]);
% ylim([p(1), p(end)]);

phonon_lost = nbar(end) - nbar(1)
figure(2)
skip = 50;                                       % Make animation faster by increasing this!
for i=1:skip:length(states)                      % Play animation  
  
%   nbar(i) = states(i).occupation_number();       % Calculate the mean occupancy
  W = states(i).wigner(X, P);                    % Calculate the wigner function
  
  subplot(2,3, [1,2,3]);
  semilogy(t(1:skip:i), nbar(1:skip:i), 'k', 'Linewidth', 2);
  hold on
  yline(nbar_env, 'r--', 'linewidth', 2);
  if nbar_0 ~= 0
    yline(nbar_0  , 'b--', 'linewidth', 2);
  end
  xlim([t(1), t(end)]);
  ylim([min([0,nbar_0-abs(phonon_lost_rough_estimate)]), nbar_0+abs(phonon_lost_rough_estimate)])
  xlabel('t [s]');
  %ylim(sort([nbar(1), nbar_final]))
  ylabel('$\bar{n}$','Interpreter','Latex', 'rotation', 0, 'VerticalAlignment', 'middle');
  hold off
  
  subplot(2,3, 5)
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

% figure(4)
% x = initial.R(1) + 5*sqrt(variance_W)*linspace(-1, +1, 150);                 % Region to plot wigner function
% p = initial.R(2) + 5*sqrt(variance_W)*linspace(-1, +1, 150);
% [X, P] = meshgrid(x, p);
% W2 = initial.wigner(X,P);
% surf(X, P, W2)
% view(0,90)
% shading interp
% xlim([x(1), x(end)])
% ylim([p(1), p(end)])
% initial.V

% ss = simulation.steady_state;
% ss.V - eye(2)*(2*nbar_env+1);



% t_inverted = linspace(0, 0.3/omega + pi/(2*omega), 1e4);                   % Timestamps for simulation
% t_harmonic = t_inverted(end) + linspace(0, pi/(2*omega), 1e4);
% 
% t = cat(2, t_inverted, t_harmonic(2:end), t_harmonic(end) + t_inverted(2:end));

% function result = omega(t)
%     omega_0 = 2*pi*305e+3;
%     
%     result = -omega_0;%*cos(4*omega_0*t);
% end

% initial = gaussian_state("coherent", 6);         % Initial state
% initial.squeeze(1.05);
% initial.displace(3 - 2*1i);
% initial.rotate(-pi/3);

%  disp( "Omega(t) = " + omega(t(i), omega_0) );