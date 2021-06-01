syms t t_0 xi s gamma omega real

A = [[    0   ,  +omega ];                 % Drift matrix
     [ -omega ,  -gamma ]];

B(t, s) = int(A, s, t);

exp_B(t,s)   = simplify( expm( B(t,s) ) );
exp_B_T(t,s) = transpose( exp_B(t,s) );

syms nbar r real
V_0 = [[(2*nbar+1)*exp(-r) ,           0       ];
       [          0        , (2*nbar+1)*exp(+r)]];

V = simplify( exp_B(t,0)*V_0*exp_B_T(t,0) )



% latex(exp_B(t,s))








% syms omega_1 omega_2 g_1 g_2 gamma_1 gamma_2 T_env  omega_cav  kappa nu_1 nu_2  s xi t t_0 real
% 
% A = [[-kappa/2 , omega_cav]; [-omega_cav   , -kappa/2]];
% 
% A_particle = [[0, omega_1]; [-omega_1, -gamma_1]];
% A = blkdiag(A, A_particle);
% 
% A_particle = [[0, omega_2]; [-omega_2, -gamma_2]];
% A = blkdiag(A, A_particle);
% 
% A(2      , 2*1+1 ) = -2*g_1*cos(nu_1*xi);
% A(2*(1+1),   1   ) = -2*g_1*cos(nu_1*xi);
% A(2      , 2*2+1 ) = -2*g_2*cos(nu_2*xi);
% A(2*(2+1),   1   ) = -2*g_2*cos(nu_2*xi);
% 
% B = int(A, xi, s, t);  % B matrix (integral of A matrix)
% 
% p = charpoly(B);       % Characteristic polynomial of B (coefs.)
% % It should be possible to evaluate the matrix exponential of B
% % using the Cayley-Hamilton theorem


