function vals = GUI_piers()
    % GUI_piers: Prompt user for ph1, ph2, BA1, BA2 in a GUI and return values

    fig = uifigure('Name', 'Enter Parameters', 'Position', [100, 100, 300, 330]);
    labels = {'pH1', 'pH2', 'SC4', 'SC6', 'PUT', 'TYR', 'pyflav'};
    n = numel(labels);
    fields = gobjects(n,1);

    for i = 1:n
        uilabel(fig, 'Text', labels{i}, 'Position', [30, 340-40*i, 60, 22]);
        fields(i) = uieditfield(fig, 'numeric', 'Position', [100, 340-40*i, 150, 22]);
    end

    btn = uibutton(fig, 'Text', 'Run', 'Position', [100, 10, 100, 30], ...
        'ButtonPushedFcn', @(btn,event) uiresume(fig));

    uiwait(fig);

    vals = zeros(1, n);
    for i = 1:n
        vals(i) = fields(i).Value;
    end

    delete(fig);
end