close all
clear 
clc

% A: front wheel center
% G: COG 
%% Parameters
load('0 - Model\parameters.mat')

%% Arbitrary parameters
% Define reference values
pman0 = 0.37 * 1e5; % 0.37 bar to Pa
n0 = 870;            % rpm
lambda_m0 = 1;       % Stoichiometric AFR

%% calculations
% Volumetric efficiency eta_v_ref
eta_v0 = eta_vn0 + eta_vn1 * n0 + eta_vn2 * n0^2 + eta_vp1 * pman0;

% Determine beta2_ref
pRatio0 = pman0 / patm;
if pRatio0 >= critpRatio
    beta20 = sqrt((2 * k / (k - 1)) * (pRatio0^(2 / k) - pRatio0^((k + 1) / k)));
else
    beta20 = sqrt(k) * (2 / (k + 1))^((k + 1) / (2 * (k - 1)));
end


% Calculate beta1_ref and alpha_ref
beta10 = (n0*Vd*eta_v0*pman0*4*(R*Tatm)^0.5)/(120*R*Tman*Ct*pi*D^2*patm*beta20);

syms alpha real
eq = (1 - cos(alpha) / cos(alpha00)) + ...
     (2 / pi) * ((a / cos(alpha)) * sqrt(cos(alpha)^2 - a^2 * cos(alpha00)^2)) + ...
     (2 / pi) * ((cos(alpha) / cos(alpha00)) * asin(a * cos(alpha00) / cos(alpha)) - ...
     a * sqrt(1 - a^2) - asin(a)) == beta10;
alpha0 = double(vpasolve(eq, alpha, [0 pi/2])) * 180 / pi; % Reference throttle angle [deg]

% Convert alpha_ref to radians for consistency
alpha0_rad = deg2rad(alpha0);

% Evaluate beta1_alpha at alpha_ref
beta1_alpha0 = (1 - cos(alpha0_rad) / cos(alpha00)) + ...
                  (2 / pi) * ((a / cos(alpha0_rad)) * ...
                  sqrt(cos(alpha0_rad)^2 - a^2 * cos(alpha00)^2)) + ...
                  (2 / pi) * ((cos(alpha0_rad) / cos(alpha00)) * ...
                  asin(a * cos(alpha00) / cos(alpha0_rad)) - ...
                  a * sqrt(1 - a^2) - asin(a));

% Define fuel film dynamics parameters
sigma3 = 1.67;
sigma4 = 0.65;
sigma5 = 9.6e-5;
sigma6 = 0.7236;

% tau_f and chi
tauf0 = sigma3 * n0^(-sigma4);
chi0 = sigma5 * n0 + sigma6;


% Calculate reference fuel injection
mfidot0 = 1000*(n0 * Vd * eta_v0 * pman0) / (lambda_m0 * R * Tman * 120 * chi0 * AFRs);
mff0 = mfidot0*(2*chi0-1);	

% State and input vector at reference
x0= [pman0; n0; mff0; lambda_m0];
u0= [mfidot0; alpha0];

y0= [n0; lambda_m0];
% e0 = [n0 - r_n0; lambda_m0 - r_lambda0];

%% Display results
% Equilibrium state
disp('x0 = ')
disp([num2str(x0(1)/1e5) ' bar'])
disp([num2str(x0(2)) ' rpm'])
disp([num2str(x0(3)) ' ?'])
disp([num2str(x0(4)) ' -'])
disp(' ')
% Equilibrium Control
disp('u0 = ')
disp([num2str(u0(1)) ' ?'])
disp([num2str(u0(2)) ' deg'])
disp(' ')