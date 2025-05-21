function param = box_gui_IDA_piers
% Define parameter names from slider_gui_IDA_piers
param_names = {'A1C', 'A1B1', 'A2B2', 'A2C', 'A2B1', 'A2B2', ...
               'lambda1', 'lambda2', 'lambda3', 'lambda4', 'lambda5', 'lambda6', ...
               'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6'};

% Create UI figure
fig = uifigure('Name', 'Parameter Inputs', 'Position', [100, 50, 400, 600]); % Adjust height

% Initialize matrix to store parameter values
num_params = numel(param_names);
param_values = zeros(num_params, 1);

% Load previous values if they exist
if isfile('piers_values.mat')
    try
        saved_data = load('piers_values.mat', 'param_values');
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
    
    % Create numeric input box
    numeric_inputs(i) = uieditfield(panel, 'numeric', ...
        'Position', [120, y_pos, 150, 22], ...
        'Value', param_values(i));
    
    % Create label for parameter name
    uilabel(panel, 'Text', param_names{i}, 'Position', [10, y_pos-5, 100, 20]);
    
    % Callback function to update values without range checking
    numeric_inputs(i).ValueChangedFcn = @(src, event) updateValue(src, i);
end

% "Run" button
btn = uibutton(fig, 'Text', 'Run', 'Position', [150, 20, 100, 30], ...
    'ButtonPushedFcn', @(btn, event) runSimulation());

% Pause execution until user presses Run
uiwait(fig);

% Nested function to update values without range checking
function updateValue(input, index)
    param_values(index) = input.Value;
end

% Nested function for Run button
function runSimulation()
    for i = 1:num_params
        param_values(i) = numeric_inputs(i).Value;
    end
    assignin('base', 'param_values', param_values); % Push to workspace
    save('piers_values.mat', 'param_values'); % Save values to file
    uiresume(fig); % Resume execution
    delete(fig); % Close GUI
end

% Return parameter names and values
param = [transpose(string(param_names)), num2cell(param_values)];

% Function to initialize default values
function defaults = initialize_defaults()
    defaults = zeros(num_params, 1);
end
end