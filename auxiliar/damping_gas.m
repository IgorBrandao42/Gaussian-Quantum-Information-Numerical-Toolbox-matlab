function gamma = damping_gas(p_gas, m_gas, T_gas, R, M)
% Calculates the damping coefficient on a nanoparticle
% Formula valid for pressures below 10mbar from PRL 116, 243601
% 
% PARAMETERS:
%   p_gas - gas pressure          [mbar]
%   R     - particle's radius     [m]
%   M     - particle's mass       [kg]
%   m_gas - mass of gas molecules [kg]
%   T_gas - gass temperature      [K]
%
% If the units above are followed, the calculated damping rate is in Hz.
%
% Take, as an example, the values from 
%
% p_aspelmeyer = 1e-6;     % [mbar]
% m_air = 28*1.67e-27;     % [kg]
% T_aspelmeyer = 300;      % [K]
% R_aspelmeyer = 71.5e-9;  % [m]
% M_aspelmeyer = 2.83e-18; % [kg]

k_B = 1.381e-23;                         % Boltzmann's Constant [J/K]

v_gas = sqrt( 3*k_B*T_gas/m_gas );       % Gas velocity from Equipartition Theorem

gamma = 15.8*R.^2*p_gas ./ ( M*v_gas );    % Damping from the environmental gas

gamma = gamma*100;                       % Convert the result to [Hz]

end



