%Matlab code by M.C. Dulce María Guerrero Tánori. If you use this code, please cite the associated article.

%Code to calculate SCP deviations in 3 subjects under different conditions (Silence, Mozart, Brahms) 
% and in 5 different states (Eyes open, musical stimulus, symbol search task, 
% paper folding and cutting task, and Stroop effect task).
% The result is printed to a .png file at 600 dpi in a 2x3 layout style.
close all; clear; clc;

fsample = 250;

output_folder = 'C:\Users\...';
if ~exist(output_folder,'dir')
    mkdir(output_folder);
end
labels5 = {'Open eyes','Music','Symbol search','Folding and cutting paper','Stroop effect'};

%Channels in the desired order
channels_selected = {'Fp1','F3','C3','P3','F7','T3','T5','O1', ...
                     'Fz','Cz','Pz','Fp2','F4','C4','P4','F8', ...
                     'T4','T6','O2'};
channels = numel(channels_selected);

%bands
freq_high = [1,4,8,12.5,1];
freq_low  = [3.5,7.5,12,35,35];

win_sec  = 10; % windows of n seconds, in this case 10 seconds
win_samp = round(win_sec*fsample);
nwin     = 20; %number of windows
rng_seed = 123; %seed

pairs_idx = nchoosek(1:5,2);
npairs = size(pairs_idx,1);

pair_names = cell(npairs,1);
for k = 1:npairs
    i = pairs_idx(k,1);
    j = pairs_idx(k,2);
    pair_names{k} = sprintf('%s vs %s', labels5{i}, labels5{j});
end

pair_colors_fixed = [
    0.00 0.45 0.74   % 1  blue
    0.85 0.33 0.10   % 2  orange
    0.93 0.69 0.13   % 3  yellow
    0.49 0.18 0.56   % 4  purple
    0.47 0.67 0.19   % 5  green
    0.30 0.75 0.93   % 6  cyan
    0.64 0.08 0.18   % 7  burgundy
    0.25 0.25 0.25   % 8  dark grey
    1.00 0.00 1.00   % 9  magenta
    0.00 0.60 0.50   % 10 teal 
];
% SILENCE_data
folders_A = {
    'C:\users...\'
};

% MOZART_data
folders_B = {
    'C:\users...\'
};

% BRAHMS_data
folders_C = {
    'C:\users...\'
};

if numel(folders_A) ~= 5 || numel(folders_B) ~= 5 || numel(folders_C) ~= 5
    error('Cada conjunto de carpetas (A, B, C) debe tener exactamente 5 carpetas.');
end

mask = triu(true(channels),1); %upper triangle

% correlation calculation block_silence
fprintf('\n==================== SILENCIO ====================\n');
Corr_5_A = zeros(channels, channels, 5);

for iState = 1:5
    folder = folders_A{iState};
    files = dir(fullfile(folder,'*.set'));
    if isempty(files)
        error('No hay .set en %s', folder);
    end

    setfile = files(1).name;
    fprintf('Procesando SILENCIO | %s (%d/5): %s\n', labels5{iState}, iState, fullfile(folder,setfile));

    EEG = pop_loadset('filename', setfile, 'filepath', folder);

    file_labels = {EEG.chanlocs.labels};
    [tf, idx] = ismember(lower(channels_selected), lower(file_labels));
    if any(~tf)
        missing = channels_selected(~tf);
        error('Faltan canales en %s: %s', setfile, strjoin(missing,', '));
    end
    EEG2 = double(EEG.data(idx,:));

    FreqBands = 5; %broand band selected 

    [B_high,A_high] = butter(4, freq_high(FreqBands)/(0.5*fsample),'high');
    EEG_high = filtfilt(B_high,A_high,EEG2'); %removes high frequencies

    [B_low,A_low] = butter(4, freq_low(FreqBands)/(0.5*fsample),'low');
    EEG_filter = filtfilt(B_low,A_low,EEG_high); %removes low frequencies

    Ref_median = median(EEG_filter,2);
    EEG_median = EEG_filter - Ref_median; %re-reference to the median

    Tlen = size(EEG_median,1);
    if Tlen < win_samp
        error('Archivo %s: no alcanza para una ventana de %d muestras.', setfile, win_samp);
    end

    maxStart = Tlen - win_samp + 1;
    rng(rng_seed + iState); %compute how many distinct points can be the start of a complete window

    if maxStart < nwin
        error('No hay suficientes posiciones (%d) para %d ventanas en %s', maxStart, nwin, setfile);
    end

    starts = randperm(maxStart, nwin).'; %generates random positions and chooses one result to be the start of a window

    % the EEG windows are extracted and concatenated into a single matrix 
    % to subsequently calculate the correlation matrix between channels
    EEG_concat = zeros(nwin*win_samp, channels);
    for w = 1:nwin
        seg = EEG_median(starts(w):starts(w)+win_samp-1, :);
        EEG_concat((w-1)*win_samp + (1:win_samp), :) = seg;
    end

    CorrMatrix = corrcoef(EEG_concat);
    CorrMatrix(1:channels+1:end) = 1;

    Corr_5_A(:,:,iState) = CorrMatrix;
end

% correlation calculation block_Mozart
fprintf('\n==================== MOZART ====================\n');
Corr_5_B = zeros(channels, channels, 5);

for iState = 1:5
    folder = folders_B{iState};
    files = dir(fullfile(folder,'*.set'));
    if isempty(files)
        error('No hay .set en %s', folder);
    end

    setfile = files(1).name;
    fprintf('Procesando MOZART | %s (%d/5): %s\n', labels5{iState}, iState, fullfile(folder,setfile));

    EEG = pop_loadset('filename', setfile, 'filepath', folder);

    file_labels = {EEG.chanlocs.labels};
    [tf, idx] = ismember(lower(channels_selected), lower(file_labels));
    if any(~tf)
        missing = channels_selected(~tf);
        error('Faltan canales en %s: %s', setfile, strjoin(missing,', '));
    end
    EEG2 = double(EEG.data(idx,:));

    FreqBands = 5;

    [B_high,A_high] = butter(4, freq_high(FreqBands)/(0.5*fsample),'high');
    EEG_high = filtfilt(B_high,A_high,EEG2');

    [B_low,A_low] = butter(4, freq_low(FreqBands)/(0.5*fsample),'low');
    EEG_filter = filtfilt(B_low,A_low,EEG_high);

    Ref_median = median(EEG_filter,2);
    EEG_median = EEG_filter - Ref_median;

    Tlen = size(EEG_median,1);
    if Tlen < win_samp
        error('Archivo %s: no alcanza para una ventana de %d muestras.', setfile, win_samp);
    end

    maxStart = Tlen - win_samp + 1;
    rng(rng_seed + iState);

    if maxStart < nwin
        error('No hay suficientes posiciones (%d) para %d ventanas en %s', maxStart, nwin, setfile);
    end

    starts = randperm(maxStart, nwin).';

    EEG_concat = zeros(nwin*win_samp, channels);
    for w = 1:nwin
        seg = EEG_median(starts(w):starts(w)+win_samp-1, :);
        EEG_concat((w-1)*win_samp + (1:win_samp), :) = seg;
    end

    CorrMatrix = corrcoef(EEG_concat);
    CorrMatrix(1:channels+1:end) = 1;

    Corr_5_B(:,:,iState) = CorrMatrix;
end

% correlation calculation block_Brahms
fprintf('\n==================== BRAHMS ====================\n');
Corr_5_C = zeros(channels, channels, 5);

for iState = 1:5
    folder = folders_C{iState};
    files = dir(fullfile(folder,'*.set'));
    if isempty(files)
        error('No hay .set en %s', folder);
    end

    setfile = files(1).name;
    fprintf('Procesando BRAHMS | %s (%d/5): %s\n', labels5{iState}, iState, fullfile(folder,setfile));

    EEG = pop_loadset('filename', setfile, 'filepath', folder);

    file_labels = {EEG.chanlocs.labels};
    [tf, idx] = ismember(lower(channels_selected), lower(file_labels));
    if any(~tf)
        missing = channels_selected(~tf);
        error('Faltan canales en %s: %s', setfile, strjoin(missing,', '));
    end
    EEG2 = double(EEG.data(idx,:));

    FreqBands = 5;

    [B_high,A_high] = butter(4, freq_high(FreqBands)/(0.5*fsample),'high');
    EEG_high = filtfilt(B_high,A_high,EEG2');

    [B_low,A_low] = butter(4, freq_low(FreqBands)/(0.5*fsample),'low');
    EEG_filter = filtfilt(B_low,A_low,EEG_high);

    Ref_median = median(EEG_filter,2);
    EEG_median = EEG_filter - Ref_median;

    Tlen = size(EEG_median,1);
    if Tlen < win_samp
        error('Archivo %s: no alcanza para una ventana de %d muestras.', setfile, win_samp);
    end

    maxStart = Tlen - win_samp + 1;
    rng(rng_seed + iState);

    if maxStart < nwin
        error('No hay suficientes posiciones (%d) para %d ventanas en %s', maxStart, nwin, setfile);
    end

    starts = randperm(maxStart, nwin).';

    EEG_concat = zeros(nwin*win_samp, channels);
    for w = 1:nwin
        seg = EEG_median(starts(w):starts(w)+win_samp-1, :);
        EEG_concat((w-1)*win_samp + (1:win_samp), :) = seg;
    end

    CorrMatrix = corrcoef(EEG_concat);
    CorrMatrix(1:channels+1:end) = 1;

    Corr_5_C(:,:,iState) = CorrMatrix;
end

% First, an average connectivity matrix is ​​calculated across the five states: SCP.
%Then, this average is subtracted from each state to determine its deviation. 
% Next, these deviations are compared pairwise using correlation. 
% Finally, these values ​​are sorted to construct a cumulative
% distribution.(SILENCE)
SCP_A = mean(Corr_5_A,3);
D_5_A = zeros(channels, channels, 5);

for iState = 1:5
    D_5_A(:,:,iState) = Corr_5_A(:,:,iState) - SCP_A;
end

R_pairs_A = zeros(npairs,1);
for k = 1:npairs
    i = pairs_idx(k,1);
    j = pairs_idx(k,2);
    v1 = D_5_A(:,:,i);
    v2 = D_5_A(:,:,j);
    R_pairs_A(k) = corr(v1(mask), v2(mask));
end

[x_sorted_A, idx_sorted_A] = sort(R_pairs_A(:));
y_A = (1:numel(x_sorted_A)).'/numel(x_sorted_A);

% First, an average connectivity matrix is ​​calculated across the five states: SCP.
%Then, this average is subtracted from each state to determine its deviation. 
% Next, these deviations are compared pairwise using correlation. 
% Finally, these values ​​are sorted to construct a cumulative
% distribution.(MOZART)
SCP_B = mean(Corr_5_B,3);
D_5_B = zeros(channels, channels, 5);

for iState = 1:5
    D_5_B(:,:,iState) = Corr_5_B(:,:,iState) - SCP_B;
end

R_pairs_B = zeros(npairs,1);
for k = 1:npairs
    i = pairs_idx(k,1);
    j = pairs_idx(k,2);
    v1 = D_5_B(:,:,i);
    v2 = D_5_B(:,:,j);
    R_pairs_B(k) = corr(v1(mask), v2(mask));
end

[x_sorted_B, idx_sorted_B] = sort(R_pairs_B(:));
y_B = (1:numel(x_sorted_B)).'/numel(x_sorted_B);

% First, an average connectivity matrix is ​​calculated across the five states: SCP.
%Then, this average is subtracted from each state to determine its deviation. 
% Next, these deviations are compared pairwise using correlation. 
% Finally, these values ​​are sorted to construct a cumulative
% distribution.(BRAHMS)
SCP_C = mean(Corr_5_C,3);
D_5_C = zeros(channels, channels, 5);

for iState = 1:5
    D_5_C(:,:,iState) = Corr_5_C(:,:,iState) - SCP_C;
end

R_pairs_C = zeros(npairs,1);
for k = 1:npairs
    i = pairs_idx(k,1);
    j = pairs_idx(k,2);
    v1 = D_5_C(:,:,i);
    v2 = D_5_C(:,:,j);
    R_pairs_C(k) = corr(v1(mask), v2(mask));
end

[x_sorted_C, idx_sorted_C] = sort(R_pairs_C(:));
y_C = (1:numel(x_sorted_C)).'/numel(x_sorted_C);

% Saving tables
T_A = table(pair_names, R_pairs_A, 'VariableNames', {'Pair','CorrDeviationMatrices'});
T_B = table(pair_names, R_pairs_B, 'VariableNames', {'Pair','CorrDeviationMatrices'});
T_C = table(pair_names, R_pairs_C, 'VariableNames', {'Pair','CorrDeviationMatrices'});

writetable(T_A, fullfile(output_folder,'R_pairs_DesvSCP_SILENCIO.csv'));
writetable(T_B, fullfile(output_folder,'R_pairs_DesvSCP_MOZART.csv'));
writetable(T_C, fullfile(output_folder,'R_pairs_DesvSCP_BRAHMS.csv'));

save(fullfile(output_folder,'R_pairs_DesvSCP_layout_3condiciones.mat'), ...
    'Corr_5_A','Corr_5_B','Corr_5_C', ...
    'SCP_A','SCP_B','SCP_C', ...
    'D_5_A','D_5_B','D_5_C', ...
    'R_pairs_A','R_pairs_B','R_pairs_C', ...
    'x_sorted_A','x_sorted_B','x_sorted_C', ...
    'idx_sorted_A','idx_sorted_B','idx_sorted_C', ...
    'pair_names','pairs_idx','pair_colors_fixed','labels5');

% Layout 2X3 figure
figure('Color','w','Position',[100 100 1500 550])

t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Plot the empirical cumulative distribution function (CDF) of the correlations
% between deviation matrices for the Silence condition. Each point corresponds
% to a pair of states, and colors indicate the specific state pair.
% The gray step line represents the cumulative distribution of all pairwise
% values.(SILENCE)
nexttile
stairs(x_sorted_A, y_A, 'LineWidth',2.5,'Color',[0.5 0.5 0.5]); hold on

for kk = 1:numel(x_sorted_A)
    this_pair_idx = idx_sorted_A(kk);
    this_color = pair_colors_fixed(this_pair_idx,:);
    plot(x_sorted_A(kk), y_A(kk), 'o', ...
        'MarkerFaceColor', this_color, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12);
end

grid on
xlim([-1 1])
ylim([0 1])
text(-0.08,1.02,'D','Units','normalized','FontSize',15,'FontWeight','bold')
xlabel('Correlation between deviation matrices')
ylabel('Cumulative proportion')

% Plot the empirical cumulative distribution function (CDF) of the correlations
% between deviation matrices for the Silence condition. Each point corresponds
% to a pair of states, and colors indicate the specific state pair.
% The gray step line represents the cumulative distribution of all pairwise
% values. (MOZART)
nexttile
stairs(x_sorted_B, y_B, 'LineWidth',2.5,'Color',[0.5 0.5 0.5]); hold on

for kk = 1:numel(x_sorted_B)
    this_pair_idx = idx_sorted_B(kk);
    this_color = pair_colors_fixed(this_pair_idx,:);
    plot(x_sorted_B(kk), y_B(kk), 'o', ...
        'MarkerFaceColor', this_color, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12);
end

grid on
xlim([-1 1])
ylim([0 1])
text(-0.08,1.02,'E','Units','normalized','FontSize',15,'FontWeight','bold')
xlabel('Correlation between deviation matrices')
ylabel('Cumulative proportion')

% Plot the empirical cumulative distribution function (CDF) of the correlations
% between deviation matrices for the Silence condition. Each point corresponds
% to a pair of states, and colors indicate the specific state pair.
% The gray step line represents the cumulative distribution of all pairwise
% values. (BRAHMS)
nexttile
stairs(x_sorted_C, y_C, 'LineWidth',2.5,'Color',[0.5 0.5 0.5]); hold on

for kk = 1:numel(x_sorted_C)
    this_pair_idx = idx_sorted_C(kk);
    this_color = pair_colors_fixed(this_pair_idx,:);
    plot(x_sorted_C(kk), y_C(kk), 'o', ...
        'MarkerFaceColor', this_color, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12);
end

grid on
xlim([-1 1])
ylim([0 1])
text(-0.08,1.02,'F','Units','normalized','FontSize',15,'FontWeight','bold')
xlabel('Correlation between deviation matrices')
ylabel('Cumulative proportion')

% Legends
h = gobjects(npairs,1);
for k = 1:npairs
    h(k) = plot(nan, nan, 'o', ...
        'MarkerFaceColor', pair_colors_fixed(k,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12);
end

lgd = legend(h, pair_names);
lgd.Layout.Tile = 'south';
lgd.NumColumns = 2;
lgd.FontSize = 9;

% Saving png 600 dpi
exportgraphics(gcf, ...
    fullfile(output_folder,'CDF_DesvSCP_layout_3condiciones_desvalSCP.png'), ...
    'Resolution', 600);