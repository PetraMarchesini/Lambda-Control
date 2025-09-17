function bar_plot_custom(modal_impact,state_labels,varargin)
    
    % Input:
    % P      - Data matrix where each row is a set of data
    % labels - Cell array of labels for each axis

    % Example data: Each row corresponds to a modal coordinate, and each column is a state variable.
    state_labels = state_labels(1:size(modal_impact,2));
    
    if ~isempty(varargin)
        % Create stacked bar plot
        b = bar(modal_impact,varargin{1});  % Stacked bar chart
    else
        b = bar(modal_impact);
    end
    
    % Custom colors for each stack (modal coordinate)
    colors = [1 0 0;  % Light blue for Mode 1
              0 1 0;  % Brown for Mode 2
              0 0 1;
              1 0 1;
              0 1 1]; % Green for Mode 3
    
    % Apply custom colors to each stack
    for k = 1:length(b)
        b(k).FaceColor = colors(k, :);  % Set the color for each stack
    end
    
    % Customize the x-axis and legend
    set(gca, 'xticklabel', state_labels,'FontWeight','bold');
end