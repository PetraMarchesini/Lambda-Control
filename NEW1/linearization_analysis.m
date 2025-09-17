%% Linearization

% Clear workspace
clear; clc;

% Define constants
k = 1.4;                    % Air isoentropic coefficient [-]
critpRatio = (2 / (k + 1))^(k / (k - 1)); % Critical pressure ratio [-]

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
alpha0 = deg2rad(5.4); % Closed throttle angle (rad)
a = 0.14;            % Geometric constant
eta_vn0 = 0.133;     % Volumetric efficiency constant
eta_vn1 = 0.391e-3;  % Volumetric efficiency constant
eta_vn2 = -0.0636e-6;% Volumetric efficiency constant
eta_vp1 = 0.202e-5;  % Volumetric efficiency constant
Hu = 43e6;           % Gasoline heating value (J/kg)
eta_b = 0.35;        % Brake thermal efficiency
I = 0.1;             % Engine moment of inertia (kg·m^2)
C = 20;              % Lambda sensor constant
Pb = 10e3;           % Brake power (W)

% Define reference values
pman0 = 0.37 * 1e5; % 0.37 bar to Pa
n0 = 870;            % rpm
lambda0 = 1;       % Stoichiometric AFR


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
eq = (1 - cos(alpha) / cos(alpha0)) + ...
     (2 / pi) * ((a / cos(alpha)) * sqrt(cos(alpha)^2 - a^2 * cos(alpha0)^2)) + ...
     (2 / pi) * ((cos(alpha) / cos(alpha0)) * asin(a * cos(alpha0) / cos(alpha)) - ...
     a * sqrt(1 - a^2) - asin(a)) == beta10;
alpha = double(vpasolve(eq, alpha, [0 pi/2])) * 180 / pi; % Reference throttle angle [deg]

% Convert alpha_ref to radians for consistency
alpha_rad = deg2rad(alpha);

% beta_1_alpha
beta1_alpha = (1 - cos(alpha) / cos(alpha0)) + ...
                  (2 / pi) * ((a / cos(alpha)) * ...
                  sqrt(cos(alpha)^2 - a^2 * cos(alpha0)^2)) + ...
                  (2 / pi) * (((cos(alpha) / cos(alpha0)) * ...
                  asin(a * cos(alpha0) / cos(alpha)) - a * sqrt(1 - a^2) - asin(a)));

% Evaluate beta1_alpha at alpha_ref
beta1_alpha0 = (1 - cos(alpha_rad) / cos(alpha0)) + ...
                  (2 / pi) * ((a / cos(alpha_rad)) * ...
                  sqrt(cos(alpha_rad)^2 - a^2 * cos(alpha0)^2)) + ...
                  (2 / pi) * ((cos(alpha_rad) / cos(alpha0)) * ...
                  asin(a * cos(alpha0) / cos(alpha_rad)) - ...
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
mfidot0 = 1000*(n0 * Vd * eta_v0 * pman0) / (lambda0 * R * Tman * 120 * chi0 * AFRs);
mffdot0 = mfidot0*(2*chi0-1);	

% State and input vector at reference
x0 = [pman0; n0; mffdot0; lambda0]; 
u0 = [mfidot0; alpha_rad];

% Define symbolic variables
syms pman n mff lambda_m mfidot alpha real



% Air mass flow (evaluated with constant beta1_alpha_ref and beta2_ref)
m_adot = (Ct * (pi / 4) * D^2 * patm / sqrt(R * Tatm)) * beta1_alpha * beta20;

% Dynamics
dpman_dt = -(n * Vd * eta_v0 * pman / (120 * V)) + (R * Tman * m_adot / V);
dn_dt = -Pb/(n*I) + (Hu * eta_b/(n*I))*mfidot;
dmff_dt = (1 / tauf0) * (-mff + chi0 * mfidot);
dlambda_m_dt = -C*lambda_m + C;

% State derivatives
dx = [dpman_dt; dn_dt; dmff_dt; dlambda_m_dt];

% Jacobians
A = jacobian(dx, [pman, n, mff, lambda_m]);
B = jacobian(dx, [mfidot, alpha]);

% Substitute constants into Jacobians
A_ref = double(subs(A, [pman, n, mff, lambda_m, mfidot, alpha]', [x0; u0]));
B_ref = double(subs(B, [pman, n, mff, lambda_m, mfidot, alpha]', [x0; u0]));

C_ref = eye(length(A_ref)); % eye=matrice identità (4per4)
D_ref = zeros(length(A_ref));

% Display results
disp('A Matrix at Reference Point:');
disp(A_ref);

disp('B Matrix at Reference Point:');
disp(B_ref);

% Display computed reference beta1_alpha
disp('Beta1_Alpha at Reference Point:');
disp(beta1_alpha0);

%% Open loop - Eigen Values / Eigen Vectors

[eigenvecA, eigenvalA] = eig(A_ref);
disp("Open Loop eigen values:")
disp(eigenvalA);
disp("Open Loop eigen vectors:")
disp(eigenvecA);

tau = abs(diag(1./eigenvalA));

%% Reachability and Observability
Co = ctrb(A_ref,B_ref);
unco = length(A_ref) - rank(Co); % Number of uncontrollable states, if 0 -> Fully reachable

disp("Reachability matrix:")
disp(Co);
if unco == 0
	disp("Linearized System fully reachable")
else
	disp("Linearized System -NOT- fully reachable")
end

Ob = obsv(A_ref,C_ref);
unob = length(A) - rank(Ob); % Number of unobservable states, if 0 -> Fully observable
disp("Observability matrix:")
disp(Ob);
if unob == 0
	disp("Linearized System fully observable")
else
	disp("Linearized System -NOT- fully observable")
end
