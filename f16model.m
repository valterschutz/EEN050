%Robust and Nonlinear Control, EEN050
%Assignment on robust controller design
%-------------------------------------------------
%-------------------------------------------------
%Initialization
clear all;
close all;
clc

%Low fidelity F16 longitudinal model from Aircraft Control and Simulation (B.L.Stevens-F.L.Lewis) pp. 156
%The linearized dynamic of the airplane F16
A=[-0.127 -235 -32.2 -9.51 0.314;-7E-4 -0.969 0 0.908 -2E-4;0 0 0 1 0; 9E-4 -4.56 0 -1.58 0;0 0 0 0 -5];
%States
%[V(ft/s) speed, alpha(rad) angle of attack, theta(rad) pitch angle, q(rad/s) pitch rate, T(lb) engine_power]'
B=[0 -0.244;0 -0.00209; 0 0;10 -0.199; 1087 0];
%Control inputs
%[thrust (N); elevator_deflection(rad)]'
C=[0 57.3 0 0 0;0 0 0 1 0;0.0208 15.2 0 1.45 0];
%Measured outputs
%[alpha(deg); q(rad/s); normal_acceleration(ft/s^2) ]
D=[0 0;0 0;0 0.033];
%Note, elevator deflection has a direct effect on vertical/normal
%accelleration
Gn=ss(A,B,C,D);

% Actuator dynamics
act_t=tf(10,[1 5 10]); % Thrust
act_e=tf(25,[1 25]); % Elevator deflection
Ga=[act_t, 0;0, act_e];
Ga = ss(Ga);

%% A1 E1
A_eigs = eig(A);
Gn_tf = tf(Gn);
ws = logspace(-2,3);
sigma_max = zeros(1,length(ws));
for i = 1:length(sigma_max)
    sigmas = svd(evalfr(Gn_tf, ws(i)));
    sigma_max(i) = max(sigmas);
end
semilogx(ws, sigma_max)
xlabel('$\omega$','Interpreter','latex')
ylabel('$\bar{\sigma}(G_n(j\omega))$', 'Interpreter','latex')
title('Maximum singular values')
Hinf = max(sigma_max);
Hinf_facit = hinfnorm(Gn);
Qo = icare(A,B,0);
H2 = sqrt(trace(B.'*Qo*B));
%% A1 E2
