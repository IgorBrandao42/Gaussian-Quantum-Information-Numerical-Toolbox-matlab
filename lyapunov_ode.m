function dVdt_vector = lyapunov_ode(~, V_old_vector, A, D)

M = size(A, 1);                            % System dimension (N_particles + 1 cavity field)partículas  + 1 campo)

V_old = reshape(V_old_vector, [M, M]);     % Vector -> matrix

dVdt = A*V_old + V_old*transpose(A) + D;   % Calculate how much the CM derivative in this time step
 
dVdt_vector = reshape(dVdt, [M^2, 1]);     % Matrix -> vector

end



