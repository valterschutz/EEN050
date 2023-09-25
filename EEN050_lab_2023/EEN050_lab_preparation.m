warning off
format compact
clc
clear all
close all

colorOrder = get(gca,'colorOrder');
%% SETUP_SRV02_EXP14_2D_GANTRY
    %
    % Sets the necessary parameters to run the SRV02 Experiment #13: Position
    % Control of 2-DOF Robot laboratory using the "s_srv02_2d_robot" and
    % "q_srv02_2d_robot" Simulink diagrams.
    % 
    % Copyright (C) 2008 Quanser Consulting Inc.
    %
    clear all;
    qc_get_step_size =1/1000; %0.001;
    deltaT = qc_get_step_size;
    Ts = qc_get_step_size;
    nPoints = (50 * 3000) -1;
%
%% Initialization Settings(DONT CHANGE ANYTHING !!)
    EXT_GEAR_CONFIG = 'HIGH';
    ENCODER_TYPE = 'E';
    TACH_OPTION = 'YES';
    LOAD_TYPE = 'NONE';
    K_AMP = 1;
    AMP_TYPE = 'VoltPAQ';
    VMAX_DAC = 10;
    ROTPEN_OPTION = '2DGANTRY-E';
    PEND_TYPE = 'MEDIUM_12IN';
    THETA_MAX = 35 * pi/180;
    ALPHA_MAX = 15.0 * pi/180;
    CONTROL_TYPE = 'AUTO';   
    X0 = pi/180*[0, 0, 0, 0];
    [ Rm, kt, km, Kg, eta_g, Beq, Jm, Jeq_noload, eta_m, K_POT, K_TACH, K_ENC, VMAX_AMP, IMAX_AMP ] = config_srv02( EXT_GEAR_CONFIG, ENCODER_TYPE, TACH_OPTION, AMP_TYPE, LOAD_TYPE );
    [ g, mp, Lp, lp, Jp_cm, Bp, RtpnOp, RtpnOff, K_POT_PEN ] = config_sp( PEND_TYPE, ROTPEN_OPTION );
    [ Lb, Jarm, K_POT_2DP, K_ENC_2DP ] = config_2d_gantry( Jeq_noload );
    K_ENC_2DIP = [-1,1].*K_ENC_2DP;
    wcf_1 = 2 * pi * 5;
    zetaf_1 = 0.9;
    wcf_2 = wcf_1;
    zetaf_2 = zetaf_1;

    %
%%%
%%%%
%%%%%
%%%%%% DO NOT CHANGE ANYTHING ABOVE THIS AREA !! Place your code below.
%%%%
%%%
%

%% Exercise 1

% Define the uncertain parameters Mp, Lp, Jp, and Co 
% using the command "ureal"

Mp = ureal('Mp', 0.1270, 'Percentage', 50);
Lp = ureal('Lp', 0.3111, 'Percentage', 50);
Jp = ureal('Jp', 0.0012, 'Percentage', 50);
Co = ureal('Co', 0.1285, 'Percentage', 10);


Lr = 0.1270; 
theta =  0;
alpha =  0;
dtheta =  0;
dalpha =  0;
Jr = 0.0083;
Dr = 0.0690;
g = 9.810;

Adelta = [0         0         1         0;
      0         0         0         1;
      0   ((Lr*Lp^3*Mp^2*dtheta^2)/2 + Lr*g*Lp^2*Mp^2)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) - ((dalpha*dtheta*Lp^4*Mp^2)/2 + 2*Jp*dalpha*dtheta*Lp^2*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)),- ((Jp*(4*Mp*alpha*dalpha*Lp^2 + 8*Dr))/2 + (Lp^4*Mp^2*alpha*dalpha)/2)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) - (- Lr*alpha*dtheta*Lp^3*Mp^2 + Dr*Lp^2*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)),-((alpha*dtheta*Lp^4*Mp^2)/2 + 2*Jp*alpha*dtheta*Lp^2*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr));
      0   (dtheta^2*(Lp^2*Lr^2*Mp^2 + Jr*Lp^2*Mp) + 2*Lp*Lr^2*Mp^2*g + 2*Jr*Lp*Mp*g)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) - (Lp^3*Lr*Mp^2*dalpha*dtheta)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)), (2*alpha*dtheta*(Lp^2*Lr^2*Mp^2 + Jr*Lp^2*Mp))/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) - (Lr*alpha*dalpha*Lp^3*Mp^2 + 2*Dr*Lr*Lp*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)),-(Lp^3*Lr*Mp^2*alpha*dtheta)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) ];
  
Bdelta = [0; 
      0;
      (4*Co*Jp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr)) + (Co*Lp^2*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr));
      (2*Co*Lp*Lr*Mp)/(Jr*Mp*Lp^2 + Jp*(4*Mp*Lr^2 + 4*Jr))];

Cdelta = [1  0  0  0;
     0  1  0  0];
 
Ddelta = [0;0];

eig(Adelta.NominalValue)

%% Exercise 2
A = blkdiag([Adelta,Adelta]);
B = blkdiag([Bdelta,Bdelta]);
C = blkdiag([Cdelta,Cdelta]);
D = blkdiag([Ddelta,Ddelta]);
P = ss(A,B,C,D);

% wiM = udyn('wiM', [6,6]);


% Parray = ufrd(P,logspace(-2,3));
Parray = usample(P, 100);
[Ppert, info] = ucover(Parray,P.NominalValue,4);
wiM = info.W1;
WiM = 10*blkdiag([wiM,wiM]);
bodemag(WiM)
%% Exercise 3

%% Exercise 4
% Design and compute the controller, and call it Chinf
% Chinf = hinfsyn(...)




%% NB: To run the simulation, the nominal model has to be in the workspace with the variable name "Pnom"
%
%%%
%%%%
%%%%%
%%%%%% Closed loop simulation environment
%%%%
%%%
%
% NB: To run the simulation, the nominal model has to be loaded to the
% workspace with the variable name "Pnom".

[ah,bh,ch,dh] = ssdata(Chinf);


figure(1)
clf;
simTime = 10;
xinit=(pi/180)*[3 3 0 0 3 3 0 0];

try
sim('Simhinf.slx')
subplot(2,1,1)
hold on
    plot(simStates.Time,simStates.Data(:,1),'linewidth',2)
    plot(simStates.Time,simStates.Data(:,3),'linewidth',2)    
    plot(simStates.Time,simStates.Data(:,2),'--','linewidth',2)
    plot(simStates.Time,simStates.Data(:,4),'--','linewidth',2)
    legend('thetaX','alphaX','thetaY','alphaY','NthetaX','NalphaX','NthetaY','NalphaY')
    ylabel('Angle [DEG]')

subplot(2,1,2)
hold on
    plot(simVoltage.Time, simVoltage.Data(:,1),'linewidth',2)
    plot(simVoltage.Time, simVoltage.Data(:,2),'--','linewidth' ,2)
    legend('Voltage X','Voltage Y')
    ylabel('Voltage [V]')
    axis([0 simTime -11 11])
catch e
    disp('Simulation failed')
end
OCL=1;

%%
%load('labb_LQR_controller.mat')
%open('ExperimentRobustControl')
