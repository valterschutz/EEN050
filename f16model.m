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
Ga_tf = tf(Ga);
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
%% Wr
Wr = [6.25^2/(s^2+2*6.25*s+6.25^2); 0; 0];
%% Wd
Wd = [0; (0.9751*s+0.2491)/(s^2+0.885*s+0.1958)];
%% We
gdc = 400;
wc = 4.3;
ghf = 0.4;
We = [makeweight(gdc,[wc,1], ghf) 0 0];
%% Wpalpha
gdc = 2.5;
wc = 0.45;
ghf = 0.015;
Wpalpha = makeweight(gdc,[wc,1], ghf);
%% Wpalphan
gdc = 2.5;
wc = 0.7;
ghf = 0.0063;
Wpan = makeweight(gdc,[wc,1], ghf);
%% Wp
Wp = [Wpalpha, 0, 0; 0, 0, Wpan];
%% Wm1
gdc = 0.2;
wc = 26;
ghf = 3;
Wm1 = makeweight(gdc,[wc,1], ghf);
% Wm2
gdc = 0.16;
wc = 42;
ghf = 2;
Wm2 = makeweight(gdc,[wc,1], ghf);
% Wm
Wm = blkdiag(Wm1, Wm2);
%% Wu
Wu = [0 1/(35*pi/180)];
%% Wn
Wn = [0.001*180/pi; 0.001; 0.001*3.28084];
%%
P_11 = [0 0; 0 0];
P_12 = [0 0 0; 0 0 0];
P_13 = Ga_tf;
P_21 = [-We*Gn_tf*Wm           ; Wp*Gn_tf*Wm       ; Wu*Wm];
P_22 = [We*Wr -We*Gn_tf*Wd 0   ; [0;0] Wp*Gn_tf*Wd [0;0]   ; 0 Wu*Wd 0];
P_23 = [-We*Gn_tf*Ga_tf        ; Wp*Gn_tf*Ga_tf    ; Wu*Ga_tf    ];
P_31 = [0 0                    ; Gn_tf*Wm         ];
P_32 = [1 0 0                  ; [0; 0; 0] Gn_tf*Wd  Wn];
P_33 = [0 0                    ; Gn_tf*Ga_tf      ];

P = [P_11, P_12, P_13;P_21 P_22 P_23;P_31 P_32 P_33];
nmeas=4;
ncont=2;
%% A2 E2
[K,CL,gamma,info] = hinfsyn(minreal(P), nmeas, ncont);
gamma
bodemag(K)
%% A2 E2
[K,CL,gamma,info] = h2syn(minreal(P), nmeas, ncont);
gamma
bodemag(K)