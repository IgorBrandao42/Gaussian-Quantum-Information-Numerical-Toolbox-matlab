A = [[  -0.1  , +2.5];
     [-2.5 ,   -0.1 ]];                             % Define the dynamics
D = zeros(2,2);
N = zeros(2,1);
t = linspace(0, 10, 500);                        % How long to simulate

initial = gaussian_state("coherent", 6);         % Initial state
initial.squeeze(1.05);
initial.displace(3 - 2*1i);
initial.rotate(-pi/3);


simulation = time_evolution(A, D, N, initial);   % Simulate!
states = simulation.run(t);

x = 1.5*linspace(-10, +10, 150);                 % Region to plot wigner function
p = 1.5*linspace(-10, +10, 150);
[X, P] = meshgrid(x, p);

for i=1:1e3                                      % Play animation  
  
  W = states(i).wigner(X, P);                    % Calculate the wigner function
  
  surf(X,P,W)                                    % Make pretty plot
  shading interp
  view(0,90)
  xlabel('x')
  xlim([x(1), x(end)])
  ylim([p(1), p(end)])
  ylabel('p')
  
  drawnow                                        % Display it
end





