%Matlab code by M.C. Dulce María Guerrero Tánori. If you use this code, please cite the associated article.

%Code to calculate SCP deviations in 3 subjects under different conditions (Silence, Mozart, Brahms) 
% and in 5 different states (Eyes open, musical stimulus, symbol search task, 
% paper folding and cutting task, and Stroop effect task).

close all; clear; clc;

fsample = 250;

output_folder = '';
if ~exist(output_folder,'dir')
    mkdir(output_folder);
end

folders_5 = {   
    '\'
    '\'
    '\'
    '\'
    '\'
};
nSubjects = numel(folders_5);
if nSubjects ~= 5
    error('folders_5 debe tener exactamente 5 carpetas.');
end

channels_selected = {'Fp1','F3','C3','P3','F7','T3','T5','O1', ...
                     'Fz','Cz','Pz','Fp2','F4','C4','P4','F8', ...
                     'T4','T6','O2'};
channels = numel(channels_selected);
channel_labels = channels_selected;

% frequency bands
Bands = {'DL','TH','AL','BT','BB'};
freq_high = [1,4,8,12.5,1];
freq_low  = [3.5,7.5,12,35,35];

win_sec  = 10; %windows of n seconds
win_samp = round(win_sec*fsample);
nwin     = 20; %number of windows
rng_seed = 123; %seed

% saving matrices
Corr_5 = zeros(channels, channels, 5);

for iSub = 1:5
    folder = folders_5{iSub};
    files = dir(fullfile(folder,'*.set'));
    if isempty(files)
        error('No hay .set en %s', folder);
    end

    setfile = files(1).name; 
    fprintf('Procesando carpeta %d: %s\n', iSub, fullfile(folder,setfile));

    EEG = pop_loadset('filename', setfile, 'filepath', folder);

    % changing order of channels
    file_labels = {EEG.chanlocs.labels};
    [tf, idx] = ismember(lower(channels_selected), lower(file_labels));
    if any(~tf)
        missing = channels_selected(~tf);
        error('Faltan canales en %s: %s', setfile, strjoin(missing,', '));
    end
    EEG2 = double(EEG.data(idx,:)); % [channels x time]

    for FreqBands = 5:5  % 
        [B_high,A_high] = butter(4, freq_high(FreqBands)/(0.5*fsample),'high');
        EEG_high = filtfilt(B_high,A_high,EEG2'); % [tiempo x canales]

        [B_low,A_low] = butter(4, freq_low(FreqBands)/(0.5*fsample),'low');
        EEG_filter = filtfilt(B_low,A_low,EEG_high); % [time x canales]

        Ref_median = median(EEG_filter,2);
        EEG_median = EEG_filter - Ref_median;

        T = size(EEG_median,1);
        if T < win_samp
            error('Archivo %s: no alcanza para una ventana de %d muestras (tiene %d).', ...
                  setfile, win_samp, T);
        end
        maxStart = T - win_samp + 1;

        rng(rng_seed + iSub);

        if maxStart < nwin
            error('No hay suficientes posiciones (%d) para %d ventanas en %s', maxStart, nwin, setfile);
        end
        starts = randperm(maxStart, nwin).';

        EEG_concat = zeros(nwin*win_samp, channels);
        for w = 1:nwin
            seg = EEG_median(starts(w):starts(w)+win_samp-1, :);  % [win_samp x channels]
            EEG_concat((w-1)*win_samp + (1:win_samp), :) = seg;   
        end

        CorrMatrix = corrcoef(EEG_concat);
        CorrMatrix(1:channels+1:end) = 1;
    end

    Corr_5(:,:,iSub) = CorrMatrix;

    writematrix(CorrMatrix, fullfile(output_folder, sprintf('CORR_%d_AL.txt', iSub)), 'Delimiter','tab');
end

Mean5 = mean(Corr_5, 3);
writematrix(Mean5, fullfile(output_folder, 'Mean_of_5_AL.txt'), 'Delimiter','tab');
save(fullfile(output_folder,'pack_5_matrices.mat'),'Corr_5','Mean5','channels_selected','win_sec','nwin','fsample','rng_seed');

% layout figure
figure('Color','w','Position',[50 50 1300 750]);
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
clim_vals = [-1 1];   

titles = { ...
    'A', 'B', 'C', ...
    'D', 'E', 'F'};

for p = 1:6
    ax = nexttile;
    
    if p <= 5
        M = Corr_5(:,:,p);
        imagesc(ax, M);
    else
        imagesc(ax, Mean5);
    end

    text(ax,-0.08,1.02,titles{p}, ...
        'Units','normalized', ...
        'FontSize',15, ...
        'FontWeight','bold');

    axis(ax,'square');
    colormap(ax, redblue);
    clim(ax, clim_vals);

    set(ax,'XTick',1:channels,'XTickLabel',channel_labels,'XTickLabelRotation',45);
    set(ax,'YTick',1:channels,'YTickLabel',channel_labels);
    set(ax,'FontSize',10);


    % cuadricula
    hold(ax,'on');
    for k = 0.5:1:channels+0.5
        plot(ax,[k k],[0.5 channels+0.5],'k-','LineWidth',0.3); %plot([k k],[0.5 channels+0.5],'k-','LineWidth',0.3);
        plot([0.5 channels+0.5],[k k],'k-','LineWidth',0.3);
    end

    % cruz Fz..Pz
    z_start = find(strcmp(channel_labels,'Fz'));
    z_end   = find(strcmp(channel_labels,'Pz'));
    if ~isempty(z_start) && ~isempty(z_end)
        plot([0.5, channels+0.5], [z_start - 0.5, z_start - 0.5], 'k-', 'LineWidth', 2);
        plot([0.5, channels+0.5], [z_end + 0.5,   z_end + 0.5],   'k-', 'LineWidth', 2);
        plot([z_start - 0.5, z_start - 0.5], [0.5, channels+0.5], 'k-', 'LineWidth', 2);
        plot([z_end + 0.5,   z_end + 0.5],   [0.5, channels+0.5], 'k-', 'LineWidth', 2);
    end
    hold(ax,'off');
end


cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Pearson correlation coefficient';
cb.Limits = clim_vals;
cb.FontSize = 14;
exportgraphics(gcf, fullfile(output_folder,'Layout_2x3_5plusMean_mozart.png'),'Resolution',600);