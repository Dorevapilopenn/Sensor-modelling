function param = slider_gui_mixed_U
    % Define parameter names
    param_names = {'H1D','H1G1','D1G1','H1G2','D1G2', ...
                   'H2D','H2G1','D2G1','H2G2','D2G2',...
                   'lambdaD', 'lambda_DH1', 'lambda_DG1','lambda_DH2', 'lambda_DG2',...
                   'H_D', 'H_DH1', 'H_DG1', 'H_DH2', 'H_DG2',...
                   'W_D', 'W_DH1', 'W_DG1','W_DH2', 'W_DG2'};

    % Create UI figure
    fig = uifigure('Name', 'Parameter Sliders', 'Position', [100, 50, 400, 1000]);

    % Define ranges
    range_0_10 = [0, 10];
    range_400_700 = [400, 700];
    range_40k_100k = [40000, 100000];
    range_50_100 = [50, 100];

    % Initialize matrix to store slider values
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

    % Panel for sliders
    panel = uipanel(fig, 'Title', 'Adjust Parameters', ...
        'Position', [10, 60, 380, 920]);

    % Store slider handles
    sliders = gobjects(num_params, 1);

    % Loop to create sliders dynamically
    for i = 1:num_params
        y_pos = 880 - (i-1) * 35;

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

        % Create slider
        sliders(i) = uislider(panel, 'Position', [120, y_pos, 150, 3], ...
            'Limits', range, 'Value', param_values(i));

        % Create label for parameter name
        uilabel(panel, 'Text', param_names{i}, 'Position', [10, y_pos-10, 100, 20]);

        % Create dynamic label to display slider value
        val_label = uilabel(panel, 'Text', num2str(sliders(i).Value, '%.2f'), ...
            'Position', [280, y_pos-10, 50, 20]);

        % Callback function to update matrix and display value
        sliders(i).ValueChangedFcn = @(src, event) updateValue(src, val_label, i);
    end

    % "Run" button
    btn = uibutton(fig, 'Text', 'Run', 'Position', [150, 20, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) runSimulation());

    % Pause execution until user presses Run
    uiwait(fig);

    % Nested function to update values
    function updateValue(slider, label, index)
        param_values(index) = slider.Value;
        label.Text = num2str(slider.Value, '%.2f');
    end

    % Nested function for Run button
    function runSimulation()
        for i = 1:num_params
            param_values(i) = sliders(i).Value;
        end
        assignin('base', 'param_values', param_values); % Optional: push to workspace
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
