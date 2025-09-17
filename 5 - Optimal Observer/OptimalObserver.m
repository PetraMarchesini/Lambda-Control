clear

addpath("0 - model\")
addpath("3 - Modes, Reachability, Observability, Open Loop\")
load("parameters3.mat")

n = size(A,1); % State dimension
% Dual system

%% Design matrices
alphad = 0.1; %1.5

% Qd matrix
s_ws = 5;  % [m/s] wind speed variance 5
s_dnum = 0.003;    % 0.04 * pi / 180; % [rad/s] gyroscope bias 0.4
s_dmu= 0.2; % [-] friction coefficient variance 0.2

Qd11 = diag([s_ws^2 s_dnum^2 s_dmu^2]);
% [s_ws^2, 0;
%         0, s_dg^2;
%         0,]; % Disturbance covariance

ld = size(Qd11,1);

s_nug = 0.0017;    % [rad/s] gyroscope noise  5 /180*pi
s_nux = 0.02;               % [m/s^2] imu noise in x 1
s_nuz = 0.02;               % [m/s^2] imu noise in z 1
s_nut = 0.07;      % [-] tonewheel noise 0.7  /60*2*pi

Qd22 = [s_nug^2, 0, 0, 0;
        0, s_nux^2, 0, 0;
        0, 0, s_nuz^2, 0;
        0, 0, 0, s_nut^2]; % Noise covariance
lnu = size(Qd22,1); 

Qd12 = zeros(ld, lnu); % Disturbance-noise covariance (correlation)

lr= size(B2,2)-(ld+lnu);

Qd = [Qd11, Qd12, zeros(ld, lr);
      Qd12.', Qd22, zeros(lnu, lr);
      zeros(lr, ld), zeros(lr, lnu), zeros(lr,lr)];

% Rd matrix (how much you don't trust sensors)
Rd = [0.00^2, 0, 0, 0;       % gyroscope 0.65
      0, 3^2, 0, 0;     % imu x 10 
      0, 0, 7^2, 0;     % imu z 10
      0, 0, 0, 0^2];      % tonewheel 0.001

Rd_bar= D2 * Qd * D2' + Rd;

%% Checks for LQR application
if isequal(Qd,Qd') && all(eig(Qd)>=0)
    disp("Qd Symmetric and Positive Semi-Definite")
else
    disp("Qd not ok")
end
if isequal(Rd,Rd') && all(eig(Rd)>=0)
    disp("Rd Symmetric and Positive Semi-Definite")
else
    disp("Rd not ok")
end
if all(eig(Rd_bar)>=0)
    disp("Rd-bar Positive Definite")
   else
    disp("Rd-bar not ok")
end
if isstabilizable((A' + alphad*eye(n)), C')
    disp("((A^T + alphad*I), C^T) Stabilizable")
   else
    disp("((A^T + alphad*I, B) not ok")
end
if isdetectable((A' + alphad*eye(n)), sqrt(Qd)*B2')
    disp("((A^T + alphad*I), sqrt(Qd)*B2^T) Detectable")
   else
    disp("((A^T + alphad*I), sqrt(Qd)*B2^T) not ok")
end
disp(' ')

%% Algebraic Riccati Equation solution
% Translating in MatLab notation
Am = A' + alphad * eye(n);
Bm = C';
Qm = B2 * Qd * B2';
Rm = Rd_bar;
Sm = (D2 * Qd * B2')';
Em = eye(n);
Gm = zeros(n);

[Xm,Km,~,info]= icare(Am,Bm,Qm,Rm,Sm,Em,Gm);
switch info.Report
    case 0
        disp("Accurate ARE solution")
    case 1
        disp("Poor accuracy ARE solution")
    case 2
        disp("Not finite ARE solution")
    case 3
        disp("No ARE solution")
end

% Optimal Observer
Ko = Km.';

% Result gain
disp('Ko = ') 
disp(num2str(Ko))