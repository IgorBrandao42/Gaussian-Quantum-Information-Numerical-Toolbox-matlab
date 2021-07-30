%% Parameters
omega = 2*pi*305e+3;                             % Particle natural frequency [Hz]
gamma = 2*pi*6.2998e-4;                          % Damping constant [Hz]

nbar_env = 3.1731e+07;                           % Environmental    occupation number
nbar_0   = 3.1731e+03;                           % Initial particle occupation number

%% Matrix definning the dynamics
A =[[    0   ,  +omega ];                        % Drift matrix for harmonic potential with damping
    [ -omega ,  -gamma ]];
        
D = diag([0, 2*gamma*(2*nbar_env+1)]);           % Diffusion matrix
N = zeros(2,1);                                  % Mean noise vector


%% Initial state
initial = gaussian_state("thermal", nbar_0);     % Thermal state with occupation number nbar_0
initial.displace(100*1i - 250);                    % Apply a displacement operator
% initial.squeeze(1.1);                            % Apply a squeezing    operator
% initial.rotate(-pi/4);                           % Apply a rotation     operator

nbar_0 = initial.occupation_number();            % Update the initial occupation number after the action of these operators


%% Timestamps for the dynamics
t = linspace(0, 5*2*pi/omega, 1e4);                % Timestamps for simulation


%% Simulation
simulation = gaussian_dynamics(A, D, N, initial);   % Simulate!
states = simulation.run(t);


%% Plot
x = 800*linspace(-1, +1, 150);                    % Region to plot wigner function
p = 800*linspace(-1, +1, 150);
[X, P] = meshgrid(x, p);

figure(2)
skip = 50;                                       % Make animation faster by skipping some timestamps!
for i=1:skip:length(states)                      % Play animation  
  
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
