function [omega, g, gamma, Delta, kappa, Gamma_recoil, Gamma_gas, E_d] = parameters(R, P_t, W_0, T_gas)

% Universal constants
c    = 2.9979*1e8;             % Speed of light             [m/s]
epsilon_0 = 8.854187817*1e-12; % Vacuum permitivity         [F/m]
h    = 6.62607004*1e-34;       % Planck constant            [J*s]
hbar = h/(2*pi);               % Reduced Planck constant    [J*s]
k_B  = 1.381e-23;                                           % Boltzmann's Constant [J/K]

% Particle
n_r = 1.475;                             % Particle refraction index
%R = 71.5*1e-9;                 % Particle radius            [m]
m_0 = 2.83*1e-18;                        % Aspelmeyer's particle mass [kg]
m = m_0*((R/(71.5*1e-9)).^3);            % ACTUAL particle mass       [kg]
alpha = 4*pi*R.^3*epsilon_0*((n_r^2-1)/(n_r^2+2)); % Particle polarizability time its volume [?]

% Tweezer
W_x = 0.67*1e-6;                         % Tweezer x waist            [m]
W_y = 0.77*1e-6;                         % Tweezer y waist            [m]
%W_0 = 0.678e-6;                          % Tweezer waist (for our calculations with gaussian beam !)
%P_t = 400*1e-3;                         % Tweezer power              [W]
e_tw = sqrt((4*P_t)/(W_x*W_y*pi*epsilon_0*c));% Tweezer amplitude [?]
k_t = 2*pi/(1064*1e-9);                  % Tweezer wavenumber         [m]
omega_tw = 2*pi*c/(1064*1e-9) ;% Tweezer frequency          [Hz]

% Cavity
L     = 1.07*1e-2;                       % Cavity length              [m]
Delta = 2*pi*315*1e3;                    % Tweezer-cavity detunnig    [Hz]
omega_cav = omega_tw + Delta;            % Cavity field frequency     [Hz]
w_0   = 41.1*1e-6;                       % Cavity waist               [m]
kappa = 2*pi*193e+3;                     % Cavity linewidth                    [Hz]
e_cav = sqrt((hbar*omega_cav)/(2*epsilon_0*pi*(w_0^2)*L)); % Cavity field amplitude [?]

% Optomechanical constants
G = (alpha*e_tw*e_cav)/(2*hbar);
G_x = (omega_cav/c)*G;

omega = sqrt( 4*alpha*P_t./(m*W_0^4*pi*epsilon_0*c) ); % Trapping frequency [Hz]
g = sqrt(hbar./(2*m.*omega)).*G_x;       % CS coupling [Hz]

% Recoil heating
I_0 = 2*P_t/(pi*W_0^2);
sigma_scatt = abs(alpha).^2*k_t^4 / (6*pi*epsilon_0^2);
P_scatt = I_0*sigma_scatt;

Gamma_recoil = P_scatt*omega_tw ./ (5*m*c^2.*omega); %[Hz]

% Thermal heating
p_gas = 1e-6;                            % Gas pressure       [mbar]
%T_gas = 130;            % Gas temperature    [K]
m_gas = 28*1.67e-27;                     % Gass molecule mass [kg]

gamma = damping_gas(p_gas, m_gas, T_gas, R, m);  % Damping [Hz]
nbar_gas = 1./( exp(hbar*omega./(k_B*T_gas)) - 1 ); % Occupation number for the gas

Gamma_gas = gamma.*nbar_gas;

E_d = sqrt( (P_t*kappa)/(hbar*omega_tw) ); % [Hz] Approximation on the cavity drive amplitude

end

% g_delic_beautiful = g_delic/(2*pi*1e+3);
% omega_beautiful   = omega  /(2*pi*1e+3);
% omega_aspelmeyer = sqrt( 4*alpha*P_t/(m*W_x^3*W_y*pi*epsilon_0*c) ); % should be: 2*pi*305*1e3 [Hz]
% W_00 = nthroot( (4*alpha*P_t)/((2*pi*305*1e3)^2*(m*pi*epsilon_0*c)) , 4); % Code to find W_0 that achieves the same natural frequency of Aspelmeyer's ground state

