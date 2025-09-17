syms t % time
syms g % Gravity Acceleration


%% Parameters

% Air
syms k critpRatio Tman R real

% Fuel
syms AFRs Hu real

% Engine
syms Vd eta_b I real

% Throttle
syms V Ct D alpha00 a real

% Volumetric Efficeincy
syms eta_vn0 eta_vn1 eta_vn2 eta_vp1 real

% Fuel Film
syms sigma3 sigma4 sigma5 sigma6 real

% Lambda sensor
syms c A om_ce real



%% System

% State
syms pman n mff lambda_m real
x= [pman; n; mff; lambda_m];

% Control
syms mfidot alpha real
u= [mfidot; alpha];

% Exogenous
syms Pb patm Tatm nu_lambda nu_n r_lambda r_n real
w= [Pb; patm; Tatm; nu_lambda; nu_n; r_lambda; r_n];

% Output
syms y  real
%y= [];

% Error
syms e  real
% e= []


%% Auxiliary equations

% Volumetric Efficiency
eta_v= eta_vn0 + eta_vn1*(n*60) +eta_vn2*(n*60)^2 +eta_vp1* pman;

% Throttle Air Massflow
beta1= (1 - cos(alpha) / cos(alpha00)) + ...
       (2 / pi) * ((a / cos(alpha)) * ...
       sqrt(cos(alpha)^2 - a^2 * cos(alpha00)^2)) + ...
       (2 / pi) * (((cos(alpha) / cos(alpha00)) * ...
       asin(a * cos(alpha00) / cos(alpha)) - a * sqrt(1 - a^2) - asin(a)));
pRatio= pman/patm;
beta2= piecewise(...
    pRatio>= critpRatio, sqrt( (2*k/(k-1))*( pRatio^(2/k) -pRatio^((k + 1)/k)) ),...
    pRatio< critpRatio, sqrt(2*k/(k+1)) * (2 / (k + 1))^(1/(k - 1))...
    );
mthdot= Ct*pi/4*D^2*patm/sqrt(R*Tatm)*beta1*beta2;

% Cylinder Air Flow
macyldot= n/2*Vd*eta_v*pman/(R*Tman);

% Fuel Film dynamics
tauf = sigma3 * (n*60)^(-sigma4);
chi = sigma5 * (n*60) + sigma6;

% Cylinder Fuel Massflow
mfdot= (1-chi)*mfidot +mff/tauf;

% Lambda
lambda= macyldot/(mfdot*AFRs);


%% Model
pman_d= -(R*Tman*macyldot)/V +(R*Tman*mthdot)/V;
n_d= -Pb/(n*4*pi^2*I) +(Hu*eta_b*mfdot)/(n*4*pi^2*I);
mffdot_d= chi*mfidot -mff/tauf;
lambda_m_d= -c*lambda_m +c*lambda;

% State Evolution xdot= f(x,u,z)
f= [pman_d;...
    n_d;...
    mffdot_d;...
    lambda_m_d];

% Output Function y= h(x,u,z)
h= [lambda_m + nu_lambda;...
    n + nu_n];

% Error Function e= he(x,u,z)
h_e= [h(1) - (1+A*sin(om_ce*t));...
      h(2) - r_n];


%% Equilibrium

% f_steady_state = solve(f==[0 0 0 0].', [mff; lambda_m; mfidot; alpha]) % not able to solve explicitely


%% Linearization

% Equilibrium triplet
syms pman_0 n_0 mff_0 lambda_m_0
x_0= [pman_0; n_0; mff_0; lambda_m_0];

syms mfidot_0 alpha_0
u_0= [mfidot_0; alpha_0];

syms Pb_0 patm_0 Tatm_0 nu_lambda_0 nu_n_0 r_lambda_0 r_n_0
w_0= [Pb_0; patm_0; Tatm_0; nu_lambda_0; nu_n_0; r_lambda_0; r_n_0];

% Performing Jacobians
A= jacobian(f,x);
B1= jacobian(f,u);
B2= jacobian(f,w);

C= jacobian(h,x);

D1= jacobian(h,u);
D2= jacobian(h,w);

Ce= jacobian(h_e,x);
De1= jacobian(h_e,u);
De2= jacobian(h_e,w);

% Calculating them in Equilibrium Triplet
A= subs(A,[x;u;w],[x_0;u_0;w_0]);
B1= subs(B1,[x;u;w],[x_0;u_0;w_0]);
B2= subs(B2,[x;u;w],[x_0;u_0;w_0]);

C= subs(C,[x;u;w],[x_0;u_0;w_0]);
D1= subs(D1,[x;u;w],[x_0;u_0;w_0]);
D2= subs(D2,[x;u;w],[x_0;u_0;w_0]);

Ce= subs(Ce,[x;u;w],[x_0;u_0;w_0]);
De1= subs(De1,[x;u;w],[x_0;u_0;w_0]);
De2= subs(De2,[x;u;w],[x_0;u_0;w_0]);

% Creating matlab functions for matrices
% linearization_A= matlabFunction(A,'File','A_fcn');
% linearization_B1= matlabFunction(B1,'File','B1_fcn');
% linearization_B2= matlabFunction(B2,'File','B2_fcn');
% 
% linearization_C= matlabFunction(C,'File','C_fcn');
% linearization_D1= matlabFunction(D1,'File','D1_fcn');
% linearization_D2= matlabFunction(D2,'File','D2_fcn');
% 
% linearization_Ce= matlabFunction(Ce,'File','Ce_fcn');
% linearization_De1= matlabFunction(De1,'File','De1_fcn');
% linearization_De2= matlabFunction(De2,'File','De2_fcn');

% Once created all these functions are merged in 1 function called linearization.m
