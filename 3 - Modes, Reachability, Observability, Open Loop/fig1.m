% Function to plot Quantities from a Single Source:
% Inputs:
%   - x= Indipendent Variable Vector
%   - y= Dependent Variable Matrix
        % Rows= Different Coordinates
        % Columns= Evolution of the Coordinate with Independent Variable
%   - Test= Name of Test in which those relationship are obtained (ex. 'Ramp Test' or 'Steady States')
        % if Test == 'Ramp Test' sets the x axis to [RPM]
        % if Test == 'Steady States' sets the x axis to [RPM]
        % if Test == 'Modal Analysis' sets the x axis to Lagrangian Coordinates
%   - Quantities= Name of Quantities that you want to plot:
        % "Displacements" , to plot Center of Mass Displacements (x,y,z)
        % "Angular Displacements" , to plot engine rotation angles (theta_x,theta_y,theta_z)
        % "Accelerations" , to plot Center of Mass Accelerations (x_dd,y_dd,z_dd)
        % "Forces" , to plot Forces on the 4 Engine Mounts (F_A,F_B,F_C,F_D)
        % "Forcing Vector Forces" , to plot Forces of the Forcing Vector (f_x, f_y, f_z)
        % "Forcing Vector Moments" , to plot Moments of the Forcing Vector (M_x, M_y, M_z)
        % "Mode Shapes" , to plot the Mode Shapes of the engine
        % Otherwise only 1st Row of y is plotted as 'Unknown'
%   - Domain= Domain in which Quantities are represented 
        % 'Time' or 'Time Histories' set x axis to [s]
        % 'Frequency' or 'Spectra' set x axis to [Hz]
        % Otherwise x axis is set to 'Unknown'
%   - Method= way in which the Quantities are obtained (ex. 'Convolution Integral' or 'Numerical Integration' or 'Max' or 'RMS')
% Output:
%   - Plot of the Quantities with 1 Subplot for each of Coordinate


function fig1(x,y,Test,Quantities,Domain,Method)
    n_eig= length(x);

    % Y Label setting through 'Quantities' name
    if Quantities == "Mode Shapes"
        ylbl= ["\lambda_1" "\lambda_2" "\lambda_3" "\lambda_4" "\lambda_5"];
        ylbl= ylbl(1:n_eig);          
    else
        ylbl= "Unknown [-]";
    end

    % X Label setting through 'Domain' name
    if (Domain == "Time") || (Domain == "Time Histories")
        xlbl= 'Time [s]';

    elseif (Domain == "Frequency") || (Domain == "Spectra")
        xlbl= 'Frequency [Hz]';

    else
        xlbl= 'Unknown [-]';
    end

    % X Label setting through 'Test' name
    if Test == "Modal Analysis"
        xlbl= ["$\theta$" "$\dot{\theta}$" "$v$" "$\omega_r$" "$\eta$"];
        xlbl= xlbl(1:n_eig);
    end

    % Finding Maximum Values to use same Scale
    if (Quantities == "Angular Displacements") || (Quantities == "Forcing Vector Moments")
        y= y(1+3:length(ylbl)+3,:);
    else
        y= y(1:length(ylbl),:);
    end
    Ymax= max(y,[],'all');
    Ymin= min(y,[],'all');

    figure('Name',[Test,' - ',Quantities])

    MrkClr= [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.4660 0.6740 0.1880]];
    LineClr= [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1];
    for i=1:length(ylbl)
        if Quantities == "Mode Shapes"
            % Subplot
            if length(ylbl) > 4
                subplot(3,2,i)
            else
                subplot(2,2,i)
            end
            
            hold on
            % plot(y(:,i),"-o","MarkerFaceColor",MrkClr(i,:),"LineWidth",1.5)
            plot(y(:,i),"-","LineWidth",1.5,"Color",LineClr(i,:))
            for j=1:length(ylbl)
                plot(j,y(j,i),'o','MarkerFaceColor',MrkClr(j,:),'MarkerSize',10)
            end
            plot(zeros(size(x)), 'k--')

            % Axis
            axis([1 n_eig 1.2*Ymin 1.2*Ymax]) % Axis Limits
            set(gca, 'FontSize',20,'FontWeight','bold') % Axis Font
            %pos= get(gca, 'Position');
            %pos(3)= pos(3) + 0.01; pos(4)= pos(4) + 0.01;
            %set(gca , 'Position', pos)

            % Labels
            x= round(x,2,'significant');
            title(append( ylbl(i)," = " ,num2str(x(i)),' [-]' ),'Color',LineClr(i,:))
            xticks(1:1:n_eig)
            xticklabels(xlbl)
            set(gca,'TickLabelInterpreter','latex')
            hold off
        else
            % Subplot
            subplot(length(ylbl),1,i);
            plot(x,y(i,:), 'LineWidth', 1.5)
    
            % Axis
            axis([min(x) max(x) 1.2*Ymin 1.2*Ymax]) % Axis Limits
            set(gca, 'FontSize',15,'FontWeight','bold') % Axis Font
    
            % Labels
            ylabel(ylbl(i), 'FontSize',15,'Rotation',0);
            xlabel(xlbl, 'FontSize',15);
            legend(Method,'FontSize',12,'FontWeight','bold')
        end

        grid on
        grid minor
    end
    set(gcf,'rend','painters','pos',[0 0 2000 1000]) % Fig position
    sgtitle(append(Test,' ',Quantities,' ',Domain,' with ',Method), 'Fontsize', 30, 'Color','red') % Title
end