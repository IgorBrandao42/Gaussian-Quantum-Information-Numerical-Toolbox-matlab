%% Input parameters
omega =    2*pi*[190e+3;  160e+3;  180e+3];      % Natural frequency of the particles  [Hz]
g     =    2*pi*[ 42e+3; 35.3e+3; 39.8e+3];      % Coupling strength                   [Hz]
gamma =    zeros( size(omega) );                 % Damping                             [Hz]
T     =    zeros( size(omega) );% ~vacuum state  % Initial temperature of each particle             [K]
T_environment = T;                               % Temperature for the environment of each particle [K]

Delta = +omega(1);                               % Cavity field natural frequency      [Hz]   (Cavity-tweezer detuning)
kappa = 0;                                       % Cavity linewidth

% Here we only consider unitary dynamics 
% of a single particle initially in the vaccum state
% in order to observe non-trivial effects on the wigner function

%% Time interval for simulation
t = linspace(0, 20/omega(1), 1e+4);              % Time stamps for the calculations


%% Example of a simple simulation only with the first particle
simple = simulation(omega(1), g(1), gamma(1), T(1), T_environment(1), Delta, kappa); % Create a simulation variable for a system with the first particle and an optical cavity

simple.run(t, "lyapunov");                       % Run the created simulation only evolving the covariance matrix in time

V_cell = simple.V_cell;                          % Get the cell with the full covariance matrices
V_cell = single_mode_CM(V_cell, 2);              % Get the covariance matrix just for the particle

x_mean = zeros(2,1);                             % Expectation value for the quadratures of the particle


%% Create a figure and define the plotting region
h = figure;
axis tight manual                                % This ensures that getframe() returns a consistent size
filename = 'test.gif';

x = linspace(-3, +3, 1e+2);                     % points in position-axis for the grid
p = linspace(-3, +3, 1e+2);                     % points in momentum-axis for the grid
[POSITION, MOMENTUM] = meshgrid(x, p);           % Create grid on phase space

for n = 1:90:length(V_cell)                      % Loop through every calculated CM
  
  V = V_cell{n};                                 % Get the current CM
  W = zeros(length(x), length(p));               % Variable to store the wigner function evaluated on the grid
  for i=1:length(p)
    X = [x; p(i)*ones(1,length(x))];
    W(:,i) = wigner(x_mean, V, X);               % Calculate the wigner function at every point on the grid
  end
  
  surf(x, p, W)                                  % Draw the surface plot
  view(0,90)
  title("Wigner function for particle 1")        % Create a title of the plot
  shading interp                                 % Make the surface with pretty colors
  xlabel('x')                                    % Label for the x axis
  ylabel('p')                                    % Label for the y axis
  xlim([p(1), p(end)])                           % Define a contant limit on the x axis
  ylim([p(1), p(end)])                           % Define a contant limit on the y axis
  %zlim([0, 3e-3])                                % Define a contant limit on the z axis (the maximum value was precalculated for this example)
  drawnow                                        % Draw the plot now
  
% Uncomment to verify that the probabillty is indeed 1 over the considered region in phase space
% area_under_wigner = dblquad(@(a,b)interp2(POSITION,MOMENTUM,W,a,b,'*cubic'), x(1), x(end),p(1), p(end))
  
% Uncomment the following lines if you wish to save the result in a gif
% frame = getframe(h);                         % Capture the plot as a frame
% im = frame2im(frame);                        % Convert the frame into a RGB image
% [imind,cm] = rgb2ind(im,256);                % Convert RGB image to indexed image
% if n == 1                                    % Write to the GIF File
%   imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
% else
%   imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 1/24);
% end
end

