function param = slider_gui_IDA
    % Define parameter names
    param_names = {'H1D','H1G1','H1G2',...
                   'H2D','H2G1','H2G2',...
                   'lambdaD', 'lambda_DH1', 'lambda_DH2',...
                   'H_D', 'H_DH1', 'H_DH2',...
                   'W_D', 'W_DH1', 'W_DH2'};

    % Create UI figure
    fig = uifigure('Name', 'Parameter Sliders', 'Position', [100, 100, 400, 650]);
    
    % Define ranges
    range_0_10 = [0, 10];         % First 6 parameters
    range_400_700 = [400, 700];   % Lambda parameters
    range_40k_100k = [40000, 100000]; % H parameters
    range_50_100 = [50, 100];     % W parameters

    % Initialize matrix to store slider values
    num_params = numel(param_names);
    param_values = zeros(num_params, 1);

    % Panel for sliders
    panel = uipanel(fig, 'Title', 'Adjust Parameters', ...
        'Position', [10, 60, 380, 570]);

    % Store slider handles
    sliders = gobjects(num_params, 1);

    % Loop to create sliders dynamically
    for i = 1:num_params %#ok<FXUP>
        y_pos = 540 - (i-1) * 35; % Adjust position dynamically

        % Determine range based on parameter type
        if i <= 6
            range = range_0_10;
        elseif contains(param_names{i}, 'lambda')
            range = range_400_700;
        elseif i <= 12
            range = range_40k_100k;
        elseif i <= 15
            range = range_50_100;
        end

        % Create slider
        sliders(i) = uislider(panel, 'Position', [120, y_pos, 150, 3], ...
            'Limits', range, 'Value', mean(range));

        % Create label for parameter name
        uilabel(panel, 'Text', param_names{i}, 'Position', [10, y_pos-10, 100, 20]);

        % Create dynamic label to display slider value
        val_label = uilabel(panel, 'Text', num2str(sliders(i).Value, '%.2f'), ...
            'Position', [280, y_pos-10, 50, 20]);

        % Store initial value
        param_values(i) = sliders(i).Value;

        % Callback function to update matrix and display value
        sliders(i).ValueChangedFcn = @(src, event) updateValue(src, val_label, i);
    end

    % "Run" button
    btn = uibutton(fig, 'Text', 'Run', 'Position', [150, 20, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) runSimulation()); %#ok<NASGU>

    % Pause execution until user presses Run
    uiwait(fig);

    % Nested function to update values
    function updateValue(slider, label, index)
        param_values(index) = slider.Value;
        label.Text = num2str(slider.Value, '%.2f');
    end

    % Nested function for Run button
    function runSimulation()
        for i = 1:num_params %#ok<FXUP>
            param_values(i) = sliders(i).Value;
        end
        assignin('base', 'param_values', param_values); % Save to workspace (optional)
        uiresume(fig); % Resume execution
        delete(fig); % Close GUI
    end
    param = [transpose(string(param_names)) param_values]
end