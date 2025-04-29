function param = box_gui_mixed_U
    % Define parameter names
    param_names = {'H1D','H1G1','D1G1','H1G2','D1G2', ...
                   'H2D','H2G1','D2G1','H2G2','D2G2',...
                   'lambdaD', 'lambda_DH1', 'lambda_DG1','lambda_DH2', 'lambda_DG2',...
                   'H_D', 'H_DH1', 'H_DG1', 'H_DH2', 'H_DG2',...
                   'W_D', 'W_DH1', 'W_DG1','W_DH2', 'W_DG2'};

    % Create UI figure
    fig = uifigure('Name', 'Parameter Inputs', 'Position', [100, 50, 400, 600]); % Adjust height

    % Define ranges
    range_0_10 = [0, 10];
    range_400_700 = [400, 700];
    range_40k_100k = [40000, 100000];
    range_50_100 = [50, 100];

    % Initialize matrix to store parameter values
    num_params = numel(param_names);
    param_values = zeros(num_params, 1);

    % Load previous values if they exist
    if isfile('slider_values_mixed.mat')
        try
            saved_data = load('slider_values_mixed.mat', 'param_values');
            if isfield(saved_data, 'param_values') && numel(saved_data.param_values) == num_params
                param_values = saved_data.param_values;
            else
                warning('Saved parameter values are incompatible. Using defaults.');
                param_values = initialize_defaults();
            end
        catch
            warning('Failed to load previous values. Using defaults.');
            param_values = initialize_defaults();
        end
    else
        param_values = initialize_defaults();
    end

    % Scrollable panel for numeric input boxes
    panel = uipanel(fig, 'Title', 'Adjust Parameters', ...
        'Position', [10, 60, 380, 500], ... % Adjust panel size
        'Scrollable', 'on'); % Enable scrolling

    % Store numeric input handles
    numeric_inputs = gobjects(num_params, 1);

    % Loop to create numeric input boxes dynamically
    for i = 1:num_params
        y_pos = 880 - (i-1) * 35; % Position elements vertically

        % Determine range based on parameter type
        if i <= 10
            range = range_0_10;
        elseif contains(param_names{i}, 'lambda')
            range = range_400_700;
        elseif i <= 20
            range = range_40k_100k;
        elseif i <= 25
            range = range_50_100;
        end

        % Create numeric input box
        numeric_inputs(i) = uieditfield(panel, 'numeric', ...
            'Position', [120, y_pos, 150, 22], ...
            'Value', param_values(i));

        % Create label for parameter name
        uilabel(panel, 'Text', param_names{i}, 'Position', [10, y_pos-5, 100, 20]);

        % Callback function to validate and update values
        numeric_inputs(i).ValueChangedFcn = @(src, event) validateAndUpdateValue(src, i, range);
    end

    % "Run" button
    btn = uibutton(fig, 'Text', 'Run', 'Position', [150, 20, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) runSimulation());

    % Pause execution until user presses Run
    uiwait(fig);

    % Nested function to validate and update values
    function validateAndUpdateValue(input, index, range)
        if input.Value < range(1) || input.Value > range(2)
            input.Value = param_values(index); % Reset to previous value
            uialert(fig, 'Value out of range!', 'Input Error');
        else
            param_values(index) = input.Value;
        end
    end

    % Nested function for Run button
    function runSimulation()
        for i = 1:num_params
            param_values(i) = numeric_inputs(i).Value;
        end
        assignin('base', 'param_values', param_values); % Push to workspace
        save('slider_values_mixed.mat', 'param_values'); % Save values to file
        uiresume(fig); % Resume execution
        delete(fig); % Close GUI
    end

    % Return parameter names and values
    param = [transpose(string(param_names)) param_values];

    % Function to initialize default values
    function defaults = initialize_defaults()
        defaults = zeros(num_params, 1);
        for j = 1:num_params
            if j <= 10
                defaults(j) = mean(range_0_10);
            elseif contains(param_names{j}, 'lambda')
                defaults(j) = mean(range_400_700);
            elseif j <= 20
                defaults(j) = mean(range_40k_100k);
            elseif j <= 25
                defaults(j) = mean(range_50_100);
            end
        end
    end
end