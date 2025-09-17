clear
close all
addpath("4 - Optimal Control\")

no_IA= load("no_integral_action2.mat");
IA= load("integral_action2.mat");

t= no_IA.out.tout;

%% Theta Plot

theta_noIA= -180/pi* no_IA.out.x.Data(:,1);
theta_IA= -180/pi* IA.out.x.Data(:,1);

figure(1)
plot(t,theta_noIA,'LineWidth',2)
hold on
plot(t,theta_IA,'LineWidth',2)

grid minor

title('Wheelie Angle (r_{\theta}= 0 deg)') %  ws= 10 m/s ; d_{\mu}= -0.2
legend('no Integral Action','Integral action')
xlabel('Time [s]')
ylabel('$\theta$ [deg]', 'interpreter', 'Latex','FontSize',15)

%% Tau_r plot

taur_noIA= no_IA.out.u.Data(:,1);
taur_IA= IA.out.u.Data(:,1);

figure(1)
plot(t,taur_noIA,'LineWidth',2)
hold on
plot(t,taur_IA,'LineWidth',2)

grid minor

title('Wheelie maneuver Control Action')
legend('Pure Stabilizer','With Integral Action')
xlabel('Time [s]')
ylabel('$\tau_r$ [Nm]', 'interpreter', 'Latex','FontSize',15)
ylim([-300 800])

% save("4 - Optimal Control\no_integral_action2.mat",'out')
%% Tracking Error plot

e_noIA= 180/pi*no_IA.out.e.Data(:,1);
e_IA= 180/pi*IA.out.e.Data(:,1);

figure(2)
plot(t,e_noIA,'LineWidth',2)
hold on
plot(t,e_IA,'LineWidth',2)

grid minor

title('Wheelie maneuver Tracking Error')
legend('Pure Stabilizer','With Integral Action')
xlabel('Time [s]')
ylabel('$\mathbf{e}$ [deg]', 'interpreter', 'Latex','FontSize',15)


% save("4 - Optimal Control\no_integral_action.mat",'out')