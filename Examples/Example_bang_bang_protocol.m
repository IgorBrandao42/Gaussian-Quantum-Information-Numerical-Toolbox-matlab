%% Parameters
omega = 2*pi*305e+3;                             % Natural frequency for the harmonic oscillator [Hz]
gamma = 0*2*pi*6.2998e-4;                        % Damping constant [Hz]
nbar_env = 2.0495e+07;                           % Environmental occupation number


%% Define the dynamics
r = 1;                                           % How much squeezing I want (e^r in the inverted potential interval
A = @(t) bang_bang(t, omega, gamma, r);          % Drift matrix (time dependent)

D = diag([0, 2*gamma*(2*nbar_env+1)]);           % Diffusion matrix
N = zeros(2,1);                                  % Mean noise vector


%% Initial state and timestamps for its time evolution
initial = gaussian_state("vacuum");              % Initial state

T = 2*r/omega + pi/(2*omega);                     % Protocol's total time
t = linspace(0, T, 1e4);                      % Timestamps for simulation (5 protocols)


%% Simulation
simulation = time_evolution(A, D, N, initial);   % Create simulation instance
states = simulation.run(t);                      % Simulate and retrieve time evolved gaussian states


%% Plot time evolving Wigner function
x = 10*linspace(-1, +1, 150);                    % Region to plot wigner function
p = 10*linspace(-1, +1, 150);
[X, P] = meshgrid(x, p);                         % 2D grid in phase-space

figure('Name', 'Bang bang protocol')             % Figure for plotting
skip = 500;                                       % Make animation faster by skipping some timestamps

for i=1:skip:length(states)                      % Play animation  
  
  A( t(i) );
  
  W = states(i).wigner(X, P);                    % Calculate the wigner function
  
  surf(X,P,W)                                    % Make pretty plot
  shading interp
  view(0,90)
  xlabel('x')
  xlim([x(1), x(end)])
  ylim([p(1), p(end)])
  ylabel('p', 'rotation', 0, 'VerticalAlignment', 'middle');
  axis square
  drawnow                                        % Display it
end


function A = bang_bang(t, omega, gamma, r)
    t_I = r/omega;                               % Inverted dynamics timespan
    t_H = pi/(2*omega);                          % Harmonic dynamics timespan
    T = 2*t_I + t_H;                             % Protocol's total time
    
    t = mod(t, T);                               % Work around times bigger than a single protocol timespan
    
    if t < t_I
      A = [[     0    , +2*omega ];              % Drift matrix for inverted potential
           [ +2*omega ,  -gamma  ]];
    elseif t <= t_I + t_H
      A = [[   0    ,  +omega ];                 % Drift matrix for harmonic potential
           [ -omega ,  -gamma ]];
    elseif t <= T
      A = [[     0    , +2*omega ];              % Drift matrix for inverted potential
           [ +2*omega ,  -gamma  ]];
    end
end

