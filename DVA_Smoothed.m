clc; clear; close all;
load("_________.mat");

cycle_values = cell_struct.equivalent_cycle_count;  % Actual cycle values
cycle_min = min(cycle_values);  % Define cycle_min
cycle_max = max(cycle_values);  % Define cycle_max

figure('Name','DVA – All Cycles','NumberTitle','off'); 
hold on; grid on;
cmap = jet(41);  % Color gradient from blue (early) to red (late)

for i = 1:41
    Q = cell_struct.AhStep_CHA{1,i};
    V = cell_struct.qOCV_CHA{1,i};
    
    if length(Q) > 10 && length(V) > 10
        % --- Downsampling ---
        N = 30;  % Balanced: high detail, low noise
        Q = Q(1:N:end);
        V = V(1:N:end);

        % --- Ensure monotonic charge curve and numeric type ---
        Q = cummax(double(Q));  % Ensure non-decreasing
        V = double(V);

        % --- Derivatives ---
        dQ = diff(Q);
        dV = diff(V);
        Q_mid = (Q(1:end-1) + Q(2:end)) / 2;

        % --- Filtering ---
        min_dQ = 1e-5;
        valid = abs(dQ) > min_dQ & abs(dV) < 0.2;
        dQ = dQ(valid);
        dV = dV(valid);
        Q_mid = Q_mid(valid);

        % --- dV/dQ and Smoothing ---
        if length(dV) > 10
            dVdQ = dV ./ (dQ + 1e-10);  % Avoid division by 0
            window = min(90, length(dVdQ));
            dVdQ_smooth = smooth(dVdQ, window, 'lowess');  % Better peak shape
            plot(Q_mid, dVdQ_smooth, 'Color', cmap(i,:), 'LineWidth', 1.5);
        end
    end
end

% --- Labels and Aesthetics ---
xlabel('Capacity (Ah)');
ylabel('dV/dQ (V/Ah)');
title('DVA – All 41 Cycles (Smoothed & Filtered)');
ylim([-1, 1]);  % Focused Y-axis for peaks
colormap(jet(41));

% Colorbar
cb = colorbar;
clim([cycle_min cycle_max]);  % Apply cycle_min and cycle_max to the color axis
cb.Ticks = linspace(cycle_min, cycle_max, 5);  % Create 5 ticks for the colorbar
cb.TickLabels = arrayfun(@(v) sprintf('%.0f', v), cb.Ticks, 'UniformOutput', false);  % Format tick labels
cb.Label.String = 'Cycle Value';
