function result = color(x)
    % Define base colors and their light variants
    colorPairs = containers.Map('KeyType', 'double', 'ValueType', 'any');

    % Store color pairs as cell arrays {base_color, light_variant}
    colorPairs(1) = {["#0984E3", "#74B9FF"]};  % Blue
    colorPairs(2) = {["#E17055", "#FAB1A0"]};  % Orange
    colorPairs(3) = {["#00CEC9", "#81ECEC"]};  % Cyan
    colorPairs(4) = {["#E84393", "#FD79A8"]};  % Pink
    colorPairs(5) = {["#00B894", "#00B894"]};  % Green (no light variant)

    % Map old indices 11-14 to their corresponding base colors (1-4)
    if x > 10
        % Return light variant of the corresponding base color
        baseIndex = x - 10;
        if baseIndex <= 4  % Only first 4 colors have light variants
            colors = colorPairs(baseIndex);
            result = colors{1}(2);  % Get light variant
        else
            error('Invalid color index');
        end
    else
        % Return base color
        if x <= 5
            colors = colorPairs(x);
            result = colors{1}(1);  % Get base color
        else
            error('Invalid color index');
        end
    end
end

