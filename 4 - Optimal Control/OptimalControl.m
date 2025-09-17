clear

addpath("0 - model\")
addpath("3 - Modes, Reachability, Observability, Open Loop\")
load("parameters3.mat")

%% Extended State

%A= Ae; B=B1e;
A=A; B=B1;

n= size(A,1);

%% Special state matrices (eps = C_eps x + D_eps u)

C_eps= eye(n); D_eps= 10000*ones(n,1);  %D_eps= zeros(n,1);
if n==4 
    %D_eps([3 4],1)=0; 
else 
    %D_eps([3 4],1)=0;
    D_eps(5,1)=0;
end

%% Design DoFs definition
% Q matrix
theta_max= -deg2rad(70) - x0(1);%-deg2rad(70) - x0(1);    %[rad] max wheelie angle
thetad_max= -deg2rad(30) - x0(2);%-deg2rad(30) - x0(2);     %[rad/s] max wheelie angular speed
v_max= 29.4 - x0(3);          %[m/s] max ideal speed of the bike 29.4 m/s
wr_max= 29.4/r - x0(4);    %[rad/s] max wheel speed
eta_max= 0.01 - 0;       %[rad s] maximum integral error (~ 5 deg for 1 second)

if n==4
    Q= 1/4*diag([1/theta_max 1/thetad_max 1/v_max 1/wr_max]).^2;
else
    Q= 1/5*diag([1/theta_max 1/thetad_max 1/v_max 1/wr_max 1/eta_max]).^2;
end

% R matrix
tau_r_max= 754.2 - u0; %[Nm]250 754.2 maximum torque at the wheel in 2nd gear

R= 1/tau_r_max^2;

% Eigenvalues worsening
alpha= 0; %0

% R bar
R_bar= R + D_eps'*Q*D_eps;

%% Checks for LQR application
if isequal(Q,Q') && all(eig(Q)>=0)
    disp("Q Symmetric and Positive Semi-Definite")
else
    disp("Q not ok")
end
if isequal(R,R') && all(eig(R)>=0)
    disp("R Symmetric and Positive Semi-Definite")
else
    disp("R not ok")
end
if all(eig(R_bar)>=0)
    disp("R-bar Positive Definite")
   else
    disp("R-bar not ok")
end
if isstabilizable((A + alpha*eye(n)), B)
    disp("((A + alpha*I), B) Stabilizable")
   else
    disp("((A + alpha*I, B) not ok")
end
if isdetectable((A + alpha*eye(n)), sqrt(Q)*C_eps)
    disp("((A + alpha*I), sqrt(Q)*C_eps) Detectable")
   else
    disp("((A + alpha*I), sqrt(Q)*C_eps) not ok")
end
disp(' ')

%% Algebraic Riccati Equation Solution
% Translation in Matlab notation
Am= A + alpha*eye(n);
Bm= B;
Qm= C_eps'*Q*C_eps;
Rm= R_bar;
Sm= (D_eps'*Q*C_eps)';
Em= eye(n);
Gm= zeros(n);

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

% Optimal Control gain
K= -Km;

% Result gain
if n==4
    disp(['Ks = [',num2str(K),']'])
else
   disp(['Ke = [',num2str(K),']'])
   disp(['Ks = [',num2str(K(:,1:4)),']'])
   disp(['Ki = [',num2str(K(:,5)),']'])
end

if n==4
    Ks= K;
else
    Ks= K(1:4);
    Ki= K(5);
end
