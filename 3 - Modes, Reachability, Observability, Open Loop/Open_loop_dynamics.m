clear
close all

addpath '0 - Model'\
load("parameters3.mat")

%% Eigenvalue Problem
% A matrix
if input("0 -- > A   ;   1 --> Ae    :        ")
    A_ol= Ae;
else
    A_ol= A;
end

% Eigenvalue problem
[eigenvectors,eigenvalues]= eig(A_ol);

% Sorting eigenvalues in descending absolute values
lambdas= diag(eigenvalues);
[lambdas,idx]= sort(abs(lambdas),'descend');
eigenvalues= eigenvalues(idx,idx);

% Reordering eigenvectors
eigenvectors= eigenvectors(:,idx);

% Modal matrix
V= eigenvectors;

% displaying values
disp('Eigenvalues:')
disp(eigenvalues)
disp(' ')
disp('Eigenvectors:')
disp(eigenvectors)

% Normalising eigenvector matrix
V_normRow= zeros(size(V));
V_normCol= zeros(size(V));
for i= 1:size(V,1) % square matrix (both dim can be used)
    V_normRow(i,:)= abs(V(i,:))/sum(abs(V(i,:)))*100;
    V_normCol(:,i)= abs(V(:,i))/sum(abs(V(:,i)))*100;
end


%% Open loop dynamics
% Choice of Mode and Time Vector
modeNum= input('Number of the mode = ');
seconds= input('Seconds = ');
t= 0:0.0001:seconds; % time vector

% Jordan Matrix
J= inv(V)*A_ol*V;

% Initial Condition
zinitial= inv(V)*([V(:,modeNum)]);
% zinitial= inv(V)*(xinitial+[-0.7854;0;0;0.01])

% Computing dynamics of modal coordinates
z= zeros(size(A_ol,1), length(t));
for time= 1:length(t)
    z(:,time)= expm(J.*t(time))*zinitial;
end

% Reverting Change of Coordinates
x= V*z;

% plotting results

figure
plot(t,x(1:4,:),'LineWidth',2)
hold on 
if size(x,1) > 4
    plot(t,x(5,:),'LineWidth',2)
    etaLabel= '$\tilde{\eta} \, [\mathrm{rad \cdot s}]$';
else
    plot(0,0,'wo')
    etaLabel= '';
end
plot(t,x(4,:)*r,"LineStyle","--",'LineWidth',2)

legend('$\tilde{\theta} \, [\mathrm{rad}]$', ...
       '$\tilde{\dot{\theta}} \, [\mathrm{rad/s}]$', ...
       '$\tilde{v} \, [\mathrm{m/s}]$', ...
       '$\tilde{\omega_r} \, [\mathrm{rad/s}]$', ... 
       etaLabel, ...
       '$\tilde{\omega_r} \cdot r \, [\mathrm{m/s}]$', ...
       'Interpreter', 'latex', 'FontSize', 13);
xlabel('Time [s]')
ylabel('State Components Amplitude','FontWeight','bold')
title(['Mode ',num2str(modeNum),' Dynamics'])
grid minor

figure
LineClr= [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1];
for i= 1:size(V,1)
    plot(t,z(i,:),'LineWidth',2,'Color',LineClr(i,:))
    hold on
end
hold off

legend('$z_1 \, [-]$', ...
       '$z_2 \, [-]$', ...
       '$z_3 \, [-]$', ...
       '$z_4 \, [-]$', ...
       '$z_5 \, [-]$', ...
       'Interpreter', 'latex', 'FontSize', 13);
xlabel('Time [s]')
ylabel('Modal Coordinates Amplitude','FontWeight','bold')
title(['Mode ',num2str(modeNum),' Dynamics'])
grid minor


fig1(diag(eigenvalues),eigenvectors,'Modal Analysis','Mode Shapes','eigenvalues','EIG')

figure
spider_plot_custom(V_normRow,{'z_1','z_2','z_3','z_4','z_5'})
if size(x,1) > 4
    legend('\theta','d\theta/dt','v','\omega_r','\eta','Location','northeast')
else
    legend('\theta','d\theta/dt','v','\omega_r','Location','northeast')
end

title('Impact of Modal Coordinates on State Variables');


figure
modal_impact = V_normRow;  % Impact in percentages
state_labels = {'\theta', 'd\theta/dt', 'v','\omega_r','\eta'};

bar_plot_custom(modal_impact,state_labels,'stacked') % ,'stacked' for stacked bars chart
legend('z_1', 'z_2', 'z_3','z_4','z_5');
ylabel('Impact [%]','FontWeight','bold');
title('Impact of Modal Coordinates on State Variables');
