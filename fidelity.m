function F = fidelity(u_1, V_1, u_2, V_2, Omega)
% Calculates the fidelity between the two arbitrary gaussian states
% \rho_1 and \rho_2. See Phys. Rev. Lett. 115, 260501.
% 
% PARAMETERS:
%   u_1 - mean value of the quadratures for state \rho_1
%   V_2 - covariance matrix             for state \rho_1
%   u_2 - mean value of the quadratures for state \rho_2
%   V_2 - covariance matrix             for state \rho_2
% 
% The user should note that non-normalized quadratures are expected;
% They are normalized to be in accordance with the notation of Phys. Rev. Lett. 115, 260501.

if size(u_1, 2) > 1                      % Just make sure there will be no problem with the dimensions !
  u_1 = u_1.';
end
if size(u_2, 2) > 1
  u_2 = u_2.';
end

N = length(u_1)/2;

if nargin == 4                           % If the symplectic form was not provided, calculate it !
  omega = [[0, 1]; [-1, 0]];             % Auxiliar variable
  Omega = [];                            % Symplectic form matrix
  for i=1:N                              % Build the symplectic form
    Omega = blkdiag(Omega, omega);
  end                                    % The option to provide the symplectic form is to optimize
end                                      % multiple calls to this method for the same system !

u_1 = u_1/sqrt(2.0);                     % Normalize the mean value of the quadratures
u_2 = u_2/sqrt(2.0);

V_1 = V_1/2.0;                           % Normalize the covariance matrices
V_2 = V_2/2.0;

delta_u = u_2 - u_1;                     % A bunch of auxiliar variables

inv_V = inv(V_1 + V_2);

V_aux = Omega.' * inv_V * ( Omega/4 + V_2*Omega*V_1 );

identity = eye(2*N);

F_tot_4 = det( 2*( sqrt(identity + (V_aux*Omega)^(-2)/4) + identity )*V_aux );

F_0 = nthroot( real(F_tot_4) / det(V_1+V_2), 4); % We take only the real part of F_tot_4 as there can be a residual complex part from numerical integration!

F = F_0*exp( -delta_u.' * inv_V * delta_u  / 4); % Fidelity

end


% if length(u_1) ~= length(u_2) || length(u_1) ~= length(V_1) || length(u_2) ~= length(V_2)
%   error("States have unmatching dimensions! Unable to calculate the fidelity between them!")
% end

