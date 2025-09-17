close all 
clear 
% clc

addpath("2 - Linearization\")

g = 9.81;               %[m/s^2] gravitational acceleration

%% Define Parameters as a Vector
% Parameters of the simulation

% Define parameters
Tman = 308;         % Intake manifold temperature (K)
AFRs = 14.67;    % Stoichiometric air fuel ratio [-]
V = 614.6e-6;        % Intake manifold volume (m^3)
Vd = 1.275e-3;      % Engine displacement (m^3)
R = 287;             % Gas constant (J/(kg·K))
Ct = 0.83;          % Flow coefficient of throttle
D = 0.05;            % Throttle bore diameter (m)
patm = 101325;      % Ambient pressure (Pa)
Tatm = 297;         % Ambient temperature (K)
alpha00 = deg2rad(5.4); % Closed throttle angle (rad)
a = 0.14;            % Geometric constant
eta_vn0 = 0.133;     % Volumetric efficiency constant
eta_vn1 = 0.391e-3;  % Volumetric efficiency constant
eta_vn2 = -0.0636e-6;% Volumetric efficiency constant
eta_vp1 = 0.202e-5;  % Volumetric efficiency constant
Hu = 43e6;           % Gasoline heating value (J/kg)
eta_b = 0.35;        % Brake thermal efficiency
I = 0.1;             % Engine moment of inertia (kg·m^2)
c = 20;              % Lambda sensor constant
Pb = 10e3;           % Brake power (W)
k = 1.4;                    % Air isoentropic coefficient [-]
critpRatio = (2 / (k + 1))^(k / (k - 1)); % Critical pressure ratio [-]
% Define fuel film dynamics parameters
sigma3 = 1.67;
sigma4 = 0.65;
sigma5 = 9.6e-5;
sigma6 = 0.7236;

%% Equilibrium Triplet (for wheelie angle of 45 deg)
% Equilibrium State
x0 = [37000;32.390037495567310;3.874997519133703e-07;1]; % [p_man, n, m_ff, lambda] 

% Equilibrium Control Action
u0 = [0.006931290720328;0.203421444918931];                                      % [kg/s] ; [rad]         

% Equilibrium Exogenous
w0 = [Pb; patm; Tatm; 0; 0; 1; 1000]; % [W]; [Pa]; [K]; [-]; [rad/s]; [-]; [rad/s];               

% Equilibrium Output
y0= [x0(2); x0(4)]; %  [-]; [rad/s]; 

% Equilibrium Output
e0= [x0(2) - w0(6); x0(4) - w0(7)]; % [-]; [rad/s];

%% Linearization
[A,B1,B2,C,D1,D2,Ce,De1,De2]= linearization(AFRs,Ct,D,Hu,I,w0(1),R,w0(3),Tman,V,Vd,a,alpha00,u0(2),c,critpRatio,eta_b,eta_vn0,eta_vn1,eta_vn2,eta_vp1,k,x0(3),u0(1),x0(2),w0(2),x0(1),sigma3,sigma4,sigma5,sigma6);


%% Initial conditions
xinitial = [0.5 * patm; 3000; 3e-5; 14]; % [p_man, n, m_ff, lambda]  initial consitions
%%%%%%%%%% DEFINE OTHER MANEUVERS

%% Simulation
simtime= 10;
timestep= 1e-4;

%% Optimal stabilizer
% Stabilizer and Integral action gains

%%%%%%%% Formal Tuning on report %%%%%%%%
Ks= [0 0 0 0;
     0 0 0 0];
Ki = [0 0 0 0;
      0 0 0 0]; 

Ke = [Ks Ki];

%% Optimal observer
% Observer gain

%%%%%%%% Formal Tuning on report %%%%%%%%
Ko = [0 0;
      0 0;
      0 0;
      0 0]; % tuning per controllo saturato

