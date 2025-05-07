clc; clear all; close all;

% Load data
load("SAR_ANR26650M1B_A_1_3.mat");

% Parameters
Cellnum = 1;
Cycnum  = numel(cell_struct.equivalent_cycle_count);
equiv   = cell_struct.equivalent_cycle_count(:);

% Create figure
figure; hold on; box on; ax = gca;
ax.FontSize = 16;
xlabel('Capacity in Ah', 'FontSize', 18);
ylabel('dV/dQ in V/Ah',  'FontSize', 18);

% Colormap settings
colormap(jet);
cmin = 0;
cmax = 1000;
caxis([cmin cmax]);

cb = colorbar;
nticks = 11;
tickVals = linspace(cmin, cmax, nticks);
cb.Ticks = tickVals;
cb.TickLabels = round(tickVals);
cb.Label.String = 'Equivalent cycle count';
cb.Label.FontSize = 16;
cb.FontSize = 14;

% Filter/smoothing settings
smoothingMethod = 'butter';   % 'butter', 'rloess', or 'wavelet'
butterOrder     = 4;
butterFc        = 0.02;
rloessWin       = 0.05;
waveletName     = 'db8';
waveletLevel    = 3;

% === DVA Loop ===
for l = 1:Cycnum
    V = cell_struct.qOCV_CHA{1,l}(:);
    Ah = cell_struct.AhStep_CHA{1,l}(:);

    if length(V) < 10 || length(Ah) < 10
        continue;
    end

    soc = Ah/Ah(end) ;  % Normalize to full charge
    Q = soc;

    % 1) Filter out non-increasing capacity
    inc_idx = [true; diff(Q) > 0];
    Q_filt  = Q(inc_idx);
    V_filt  = V(inc_idx);

    % 2) Unique capacity values
    [Quniq, idxQ] = unique(Q_filt);
    Vuniq = V_filt(idxQ);

    % 3) Smoothing
    switch lower(smoothingMethod)
        case 'butter'
            fs = 1 / mean(diff(Quniq));
            [b, a] = butter(butterOrder, butterFc, 'low');
            V_smooth = filtfilt(b, a, Vuniq);
        case 'rloess'
            V_smooth = smoothdata(Vuniq, 'rloess', ...
                        floor(rloessWin * numel(Vuniq)));
        case 'wavelet'
            V_smooth = wdenoise(Vuniq, waveletLevel, ...
                        'Wavelet', waveletName, ...
                        'DenoisingMethod', 'SURE');
        otherwise
            error('Unknown smoothing method.');
    end

    % 4) Differentiate dV/dQ
    dVdQ_raw = diff(V_smooth) ./ diff(Quniq);

    % 5) Z-Score Outlier Removal (Thresholding)
    % Remove NaN values
    dVdQ_raw_noNaN = dVdQ_raw(~isnan(dVdQ_raw));

    % Calculate Z-score
    zscore_dVdQ = (dVdQ_raw_noNaN - mean(dVdQ_raw_noNaN)) / std(dVdQ_raw_noNaN);

    % Set a Z-score threshold to remove extreme values
    z_threshold = 3;  % You can lower this value (e.g., 2.5) to remove more outliers
    outlier_idx = abs(zscore_dVdQ) > z_threshold;  % Find the indices of the extreme peaks
    dVdQ_raw(outlier_idx) = NaN;  % Replace outliers with NaN

    % Fill missing values and ensure finite input for filtfilt
    dVdQ_filt = fillmissing(dVdQ_raw, 'movmedian', 20);  % Increase window for smoother result
    dVdQ_filt(~isfinite(dVdQ_filt)) = 0;  % Set remaining non-finite values to 0

    % Optional final smoothing (Increase the order or window size for more smoothing)
    dVdQ = filtfilt(ones(1,3000)/3000, 1, dVdQ_filt);  % Increased filter window size for stronger smoothing
    %dVdQ = -dVdQ;
    % 6) Plot
    xPlot = Quniq(1:end-1);
    yPlot = dVdQ;

    cmap = jet(256);
    cv = (equiv(l) - cmin) / (cmax - cmin);
    idx = max(1, min(256, round(cv * 255) + 1));
    clr = cmap(idx, :);

    plot(xPlot, yPlot, 'LineWidth', 2, 'Color', clr);
end

title(['DVA – LFP50 Cell0' num2str(Cellnum)], 'FontSize', 20);
set(gcf, 'Position', [100, 100, 800, 550]);
ylim([-1 1]);  % Typical for LFP cells
