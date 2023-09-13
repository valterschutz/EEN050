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
loglog(ws, sigma_max)
xlabel('$\omega$','Interpreter','latex')
ylabel('$\bar{\sigma}(G_n(j\omega))$', 'Interpreter','latex')
title('Maximum singular values')
grid()
set(gca, 'FontSize', 24);
h2norm = norm(Gn,2);
hinfnorm = norm(Gn,'Inf');
%% A1 E2
Q = diag([0,1,0,0,0,1]);
R = diag([1,1]);
Ae = [A zeros(5,1); -1 zeros(1,5)];
Be = [B; zeros(1,2)];
Ce = [C zeros(3,1)];
Gn_e = ss(Ae,Be,Ce,D);
K = lqr(Gn_e,Q,R);
Gc = (A-B*K(:,1:5));
%% A2
s = tf("s");
%% We
gdc = 400;
wc = 4.3;
ghf = 0.4;
[k,z,p] = getparams(gdc,ghf,wc);
We = k * (s + z)/(s + p)
%% Wpalpha
gdc = 2.5;
wc = 0.45;
ghf = 0.015;
[k,z,p] = getparams(gdc,ghf,wc);
Wealpha = k * (s + z)/(s + p)
%% Wpan
gdc = 2.5;
wc = 0.7;
ghf = 0.0063;
[k,z,p] = getparams(gdc,ghf,wc);
Wpan = k * (s + z)/(s + p)
%% Wm1
gdc = 0.2;
wc = 26;
ghf = 3;
[k,z,p] = getparams(gdc,ghf,wc);
Wm1 = k * (s + z)/(s + p)
%% Wm2
gdc = 0.16;
wc = 42;
ghf = 2;
[k,z,p] = getparams(gdc,ghf,wc);
Wm2 = k * (s + z)/(s + p)