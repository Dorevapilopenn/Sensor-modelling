clc; clear; close all;

%% Create GUI for User Inputs
fig = figure('Name', 'Equilibrium Simulation', 'NumberTitle', 'off', 'Position', [100 100 1000 700]);

% Parameter Names
param_names =  {'H1D1','H1G1','H1G2',...
    'H2D2','H2G1','H2G2',...
    'lambdaD', 'H_D', 'W_D',...
    'lambda_DH1', 'H_DH1', 'W_DH1',...
    'lambda_DH2', 'H_DH2', 'W_DH2'};
num_params = length(param_names);
sliders = zeros(1, num_params);  % Store sliders
text_values = zeros(1, num_params);  % Store text handles for real-time updates

% Create sliders and labels
for i = 1:num_params
    uicontrol('Style', 'text', 'Position', [50, 650-50*i, 80, 20], 'String', param_names{i});
    sliders(i) = uicontrol('Style', 'slider', 'Position', [140, 650-50*i, 200, 20], ...
        'Min', 1, 'Max', 10, 'Value', 5, 'Callback', @(src,~) update_values());
    text_values(i) = uicontrol('Style', 'text', 'Position', [350, 650-50*i, 80, 20], 'String', num2str(get(sliders(i), 'Value')));
end

% Gaussian Absorption Slider for D
uicontrol('Style', 'text', 'Position', [400, 600, 150, 20], 'String', 'D Gaussian Width');
slider_D_width = uicontrol('Style', 'slider', 'Position', [550, 600, 200, 20], ...
    'Min', 50, 'Max', 100, 'Value', 70, 'Callback', @(src,~) update_values());

uicontrol('Style', 'text', 'Position', [400, 550, 150, 20], 'String', 'D Gaussian Peak');
slider_D_peak = uicontrol('Style', 'slider', 'Position', [550, 550, 200, 20], ...
    'Min', 400, 'Max', 600, 'Value', 480, 'Callback', @(src,~) update_values());

% Button to Run Simulation
uicontrol('Style', 'pushbutton', 'Position', [350, 50, 100, 40], 'String', 'Run', 'FontSize', 12, ...
    'Callback', @(src,~) run_simulation());

% Store GUI values (now update when the button is clicked)
function update_values()
    % Update slider values and real-time text labels
    for j = 1:num_params
        param_values(j) = get(sliders(j), 'Value');
        set(text_values(j), 'String', num2str(param_values(j)));  % Update real-time display
    end
end

%% Simulation Function
function run_simulation()
    % Get slider values for the parameters
    param_values = zeros(1, num_params);
    for j = 1:num_params
        param_values(j) = get(sliders(j), 'Value');
    end

    % Randomize Gaussian parameters for other species
    nsamp = 30;  % Number of simulations
    lam = 400:1:600;  % Wavelength range
    num_species = 6;

    % Absorption parameters
    mean_vals = [get(slider_D_peak, 'Value'), randi([420 520], 1, num_species - 1)];
    width_vals = [get(slider_D_width, 'Value'), 50 + 30 * rand(1, num_species - 1)];
    height_vals = [50000, 50000 + 50000 * rand(1, num_species - 1)];

    % Generate absorption spectra
    A = zeros(num_species, length(lam));
    for i = 1:num_species
        A(i, :) = height_vals(i) * exp(-log(2) * ((lam - mean_vals(i)) ./ (width_vals(i)/2)).^2);
    end

    % Equilibrium Calculation (Newton-Raphson)
    Model = [1 0 0 1 1 0; 0 1 0 1 0 1; 0 0 1 0 1 1];
    beta = 10.^param_values;
    C_tot = [16.5e-6 * ones(nsamp, 1), 0.0001 * rand(nsamp, 1), 16.5e-6 * ones(nsamp, 1)];
    c_guess = [1e-10 1e-10 1e-10];  % Initial guess for concentrations

    % Initialize concentration matrix C for storing results
    C = zeros(nsamp, num_species);
    
    % Run Newton-Raphson method for each simulation
    for i = 1:nsamp
        C(i, :) = NewtonRaphson(Model, beta, C_tot(i, :), c_guess);
        c_guess = C(i, 1:3);  % New guess for next iteration
    end
    
    % Aggregate absorption profiles
    C_colored = C(:, [3, 5, 6]);  % Select absorbing species
    D_calc = C_colored * A;

    % Perform PCA
    [coeff, score, ~] = pca(D_calc);

    % Plot PCA Results
    subplot(1, 1, 1);
    scatter(score(:, 1), score(:, 2), 30, 'filled');
    xlabel('PC1');
    ylabel('PC2');
    title('PCA Results (30 Simulations)');
end

%% Newton-Raphson Function
function c_spec = NewtonRaphson(Model, beta, c_tot, c)
    ncomp = length(c_tot);
    nspec = length(beta);
    c_tot(c_tot == 0) = 1e-15;  % Avoid division by zero

    it = 0;
    while it <= 99
        it = it + 1;
        c_spec = beta .* prod(repmat(c', 1, nspec) .^ Model, 1);
        c_tot_calc = sum(Model .* repmat(c_spec, ncomp, 1), 2)';
        d = c_tot - c_tot_calc;
        
        if all(abs(d) < 1e-15)
            return;
        end
        
        for j = 1:ncomp
            for k = j:ncomp
                J_s(j, k) = sum(Model(j, :) .* Model(k, :) .* c_spec);
                J_s(k, j) = J_s(j, k);  % Symmetric Jacobian
            end
        end
        
        delta_c = (d / J_s) * diag(c);
        c = c + delta_c;
        
        while any(c <= 0)
            delta_c = 0.5 * delta_c;
            c = c - delta_c;
            if all(abs(delta_c) < 1e-15)
                break;
            end
        end
    end
    
    if it > 99
        fprintf(1, 'No convergence\n');
    end
end
