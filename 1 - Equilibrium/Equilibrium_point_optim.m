close all
clear 
clc

% cd ..
% addpath(cd)

%% Parameters
load('0 - Model\parameters.mat')

%% Arbitrary parameters
pman0= 0.37 * 1e5; % 0.37 bar to Pa
lambda_m0= 1; 

w0 = [Pb; patm; Tatm; 0; 0; 1; 1000];
r_n0= w0(end-1);
r_lambda0= w0(end);

%% Calculation
% ---------------------------------------- setup
% Optimization Variables
n0    = optimvar('n0',1,'LowerBound',0); %[-5;35]
mff0   = optimvar('mff0',1,'LowerBound',0);
mfidot0 = optimvar('mfidot0',1,'LowerBound',0); % 
alpha0=  optimvar('alpha0',1,'LowerBound',0,'UpperBound',deg2rad(90)); %[-5;35]

% Starting points
IC.n0    = 870/60;
IC.mff0   = 1e-5;
IC.mfidot0 = 0.1633;
IC.alpha0= alpha00;

% --------------------------------------------- Define equations
% Volumetric Efficiency
eta_v= eta_vn0 + eta_vn1*(n0*60) +eta_vn2*(n0*60)^2 +eta_vp1* pman0;

% Throttle Air Massflow
beta1= (1 - cos(alpha0) / cos(alpha00)) + ...
       (2 / pi) * ((a / cos(alpha0)) * ...
       sqrt(cos(alpha0)^2 - a^2 * cos(alpha00)^2)) + ...
       (2 / pi) * (((cos(alpha0) / cos(alpha00)) * ...
       asin(a * cos(alpha00) / cos(alpha0)) - a * sqrt(1 - a^2) - asin(a)));
pRatio= pman0/patm;
if pRatio>= critpRatio
    beta2 = sqrt( (2*k/(k-1))*( pRatio^(2/k) -pRatio^((k + 1)/k)) );
else
    beta2 = sqrt(2*k/(k+1)) * (2 / (k + 1))^(1/(k - 1));
end
mthdot= Ct*pi/4*D^2*patm/sqrt(R*Tatm)*beta1*beta2;

% Cylinder Air Flow
macyldot= n0/2*Vd*eta_v*pman0/(R*Tman);

% Fuel Film dynamics
tauf = sigma3 * (n0*60)^(-sigma4);
chi = sigma5 * (n0*60) + sigma6;

% Cylinder Fuel Massflow
mfdot= (1-chi)*mfidot0 +mff0/tauf;

% Lambda
lambda= macyldot/(mfdot*AFRs);


% Model
pman_d= -(R*Tman*macyldot)/V +(R*Tman*mthdot)/V;
n_d= -Pb/(n0*4*pi^2*I) +(Hu*eta_b*mfdot)/(n0*4*pi^2*I);
mffdot_d= chi*mfidot0 -mff0/tauf;
lambda_m_d= -c*lambda_m0 +c*lambda;

% ---------------------------------------------- Optimization Problem
% The goal is to find the values of v0, wr0, taur0 that keep the bike at
% equilibrium when the wheeling angle is equal to theta0.
prob = eqnproblem;
% equation 1
eq1 = pman_d == 0 ;
prob.Equations.eq1 = eq1;
% equation 2
eq2 = n_d == 0 ;
prob.Equations.eq2 = eq2;
% equation 3
eq3 = mffdot_d == 0 ;
prob.Equations.eq3 = eq3;
% equation 4
eq4 = lambda_m_d == 0 ;
prob.Equations.eq4 = eq4;

sol = solve(prob,IC);

% Solution
n0 = sol.n0;
mff0 = sol.mff0;
mfidot0 = sol.mfidot0;
alpha0= sol.alpha0;

x0= [pman0; n0; mff0; lambda_m0];
u0= [mfidot0; alpha0];

y0= [n0; lambda_m0];
e0 = [n0 - r_n0; lambda_m0 - r_lambda0];

%% Display results
% Equilibrium state
disp('x0 = ')
disp([num2str(x0(1)/1e5) ' bar'])
disp([num2str(x0(2)*60) ' rpm'])
disp([num2str(x0(3)*1e6) ' mg'])
disp([num2str(x0(4)) ' -'])
disp(' ')
% Equilibrium Control
disp('u0 = ')
disp([num2str(u0(1)*1e3) ' g/s'])
disp([num2str(u0(2)*180/pi) ' deg'])
disp(' ')

% rmpath(cd)