function spider_plot_custom(P, labels)
    % Input:
    % P      - Data matrix where each row is a set of data
    % labels - Cell array of labels for each axis

    % Set labels length based on data
    labels= labels(1:size(P,2));

    % Determine the number of data sets and axes
    [num_data, num_vars] = size(P);
    
    % Calculate the angle for each axis
    theta = linspace(0, 2*pi, num_vars + 1);  % One extra point to close the loop
    
    % Create a polar axes object
    polaraxes;  % Creates polar axes
    hold on
    
    % Plot each data set
    MrkClr= [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1];
    for i = 1:num_data
        data = [P(i, :), P(i, 1)];  % Repeat the first value to close the spider
        polarplot(theta, data, '-', 'LineWidth', 2)  % Polar plot for each row
    end
    for i=1:length(labels)
        data = [P(i, :), P(i, 1)];  % Repeat the first value to close the spider
        for j=1:length(labels)
            polarplot(theta(j), data(j), 'o', 'MarkerFaceColor', MrkClr(j,:),'MarkerSize',10)  % Polar plot for each row
        end
    end
    
    % Customize the plot limits
    rlim([0 max(P(:))])  % Adjust the radial limits

    % Customize the grid and labels
    ax = gca;
    ax.ThetaTick = rad2deg(theta(1:end-1));  % Set angular ticks
    ax.ThetaTickLabel = labels;  % Assign axis labels
    ax.FontWeight= "bold";
    ax.FontSize= 13;
    ax.RAxisLocation = 0;  % Set radial axis location
    ax.RTickLabel= {};
    ax.ThetaGrid = 'on';  % Show angular grid
    ax.RGrid = 'on';  % Show radial grid
    ax.LineWidth = 1.5;  % Set grid line thickness
    %ax.RAxis.Label.String = "Impact [%]";
    
    hold off
end