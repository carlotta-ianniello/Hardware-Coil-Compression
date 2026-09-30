% =========================================================
% pca_compression_pipeline.m
%
% PCA-based channel compression pipeline for a 4-channel
% hardware mixing matrix.
%
% Workflow:
%   1.  Load signal and noise data
%   2.  Compute noise covariance Psi (for pre-whitening)
%   3.  If n_coils=8: evaluate all candidate 4-channel subsets
%       and select the one whose compression matrix phases best
%       fit the hardware phase range, with priority on modes 1 & 2
%   4.  Build ROI mask (SNR threshold, radial crop, or both)
%   5.  Compute pre-whitened PCA via SVD
%   6.  Optimise per-mode global phase rotation into hardware range
%   7.  Compute SNR retention maps for 4/3/2/1 virtual modes
%   8.  Save compression matrix A for calculate_settings_complex.m
%
% Output:
%   compression_matrix.mat  ?  variable A [4 x 4] complex,
%                              ready for loadC path in
%                              calculate_settings_complex.m
% =========================================================

addpath(genpath('C:\Users\iannic01\Personal\Matlab_tools\Coil_compression'))
addpath(genpath('C:\Users\iannic01\Personal\Matlab_tools\SNR_toolbox_Ryan_1.5T'))
addpath(genpath('C:\Users\iannic01\Personal\Matlab_tools\mrir'))
addpath('C:\Users\iannic01\Personal\Matlab_tools')
clear; clc;

%% ========================
%  CONFIGURATION
%% ========================

% --- Raw data format ---
% 'vd'     : VD13+ Siemens format, read via mapVBVD
% 'legacy' : older format, read via read_meas_dat_UNIX
data_format  = 'legacy';   % 'vd' | 'legacy'
noise_source = 'file';     % 'prescan' | 'file' | 'kspace'
kspace_noise_fraction = 0.05;   % fraction of outer ky lines used when
                                  % noise_source = 'kspace'

% --- Hardware phase range (degrees, from calibration output) ---
hw_ph_lo  = -68.4;
hw_ph_hi  = +117.2;
hw_ph_buf =  10.0;

% --- PCA options ---
prewhiten       = true;
roi_type = 'manual';   % 'snr' | 'radial' | 'both' | 'manual'
snr_threshold   = 0.10;
radial_fraction = 2/3;
save_roi        = 0;
load_roi        = 1;
roi_filename    = 'roi_mask3.mat';


% --- Coil subset selection weights (n_coils = 8 only) ---
mode_weights = [4.0, 3.0, 1.0, 0.5];

% --- Manual channel override (n_coils = 8 only) ---
% Set manual_channels to a 1x4 vector to bypass automatic subset
% selection entirely, e.g. manual_channels = [3 4 5 6].
% Set to [] to use automatic selection.
manual_channels = [7 8 1 2];

% --- Output ---
save_compression_matrix = 0;
output_file  = 'compression_matrix_roi3.mat';
plot_results = true;

%% ========================
%  STEP 1: FILE SELECTION
%% ========================
fprintf('=== PCA Compression Pipeline ===\n\n');
fprintf('Step 1: Select raw data files\n');

% Signal scan
[sig_name, sig_path] = uigetfile( ...
    {'*.dat;*.out', 'Raw data (*.dat, *.out)'; '*.*', 'All Files'}, ...
    'Select SIGNAL scan (phantom / subject)');
if isequal(sig_name, 0), error('No signal file selected.'); end
signal_filepath = fullfile(sig_path, sig_name);
fprintf('  Signal : %s\n', sig_name);

% Noise file: only requested when noise_source = 'file'
noise_filepath = '';
if strcmp(noise_source, 'file')
    [nse_name, nse_path] = uigetfile( ...
        {'*.dat;*.out', 'Raw data (*.dat, *.out)'; '*.*', 'All Files'}, ...
        'Select NOISE scan (RF off)', sig_path);
    if isequal(nse_name, 0), error('No noise file selected.'); end
    noise_filepath = fullfile(nse_path, nse_name);
    fprintf('  Noise  : %s\n', nse_name);
end

%% ========================
%  STEP 2: LOAD SIGNAL DATA
%% ========================
fprintf('\nStep 2: Loading signal scan...\n');
scan_params = struct();

if strcmp(data_format, 'vd')
    [images, kspace_sig, scan_params] = load_pca_signal_vd(signal_filepath, scan_params, config);
else
    [images, kspace_sig, scan_params] = load_pca_signal_legacy(signal_filepath, scan_params);
end

% Handle multiple-acquisition averaging.
% mapVBVD may return data as [Nx, Ny, Nz, Ncoil, 1, 1, ..., Navg] when
% the sequence recorded more than one average and averaging was not applied
% in the loader.  Average coherently in image space (equivalent to k-space
% averaging since FFT is linear over complex data), then squeeze any
% singleton dimensions before standardising to 4-D.
if ndims(images) > 4
    avg_dim = ndims(images);       % averages are always the last dimension
    n_avg   = size(images, avg_dim);
    if n_avg > 1
        images = mean(images, avg_dim);
        fprintf('  Averaged %d acquisitions\n', n_avg);
    end
    images = squeeze(images);      % remove remaining singletons
end

% Enforce 4-D: [Nx, Ny, Nz, Ncoil]
if ndims(images) == 3
    images = reshape(images, size(images,1), size(images,2), 1, size(images,3));
end
[Nx, Ny, Nz, n_coils] = size(images);
fprintf('  Images : %d x %d x %d,  %d coils\n', Nx, Ny, Nz, n_coils);

if ~(n_coils == 4 || n_coils == 8 || n_coils == 9)
    error('Expected 4, 8, or 9 coils; found %d.', n_coils);
end

% 9-channel coil: channel 1 is a dummy.
% Strip it from images now so n_coils is correct before noise loading.
% noise_raw is stripped immediately after it is loaded (Step 3 below).
dummy_ch = (n_coils == 9);
if dummy_ch
    images  = images(:, :, :, 2:end);
    n_coils = 8;
    fprintf('  9-channel coil: channel 1 (dummy) removed ? %d coils\n\n', n_coils);
end

%% ========================
%  STEP 3: LOAD NOISE DATA
%% ========================
fprintf('\nStep 3: Loading noise...\n');

switch noise_source

    case 'prescan'
        % VD format: twix.noise() contains a dedicated noise pre-scan
        % recorded by the scanner before the imaging sequence.
        % This is the preferred source: it captures the actual noise
        % environment including all amplifier and coupling effects.
        if ~strcmp(data_format, 'vd')
            warning('prescan noise requires VD format. Falling back to kspace.');
            noise_source = 'kspace';
        else
            try
                twix = mapVBVD(signal_filepath);
                if ~isfield(twix, 'noise') || isempty(twix.noise)
                    warning('No noise pre-scan in this file. Falling back to kspace.');
                    noise_source = 'kspace';
                else
                    noise_raw = double(squeeze(twix.noise()));
                    % noise() returns [Nsamples x Ncoil] or [Ncoil x Nsamples]
                    if size(noise_raw, 2) ~= n_coils && size(noise_raw, 1) == n_coils
                        noise_raw = noise_raw.';   % ensure [Nsamples x Ncoil]
                    end
                    fprintf('  Source : scanner pre-scan  (%d samples x %d coils)\n', ...
                        size(noise_raw,1), size(noise_raw,2));
                end
            catch ME
                warning('Failed to read noise pre-scan: %s. Falling back to kspace.', ME.message);
                noise_source = 'kspace';
            end
        end

    case 'file'
        % Separate noise-only scan (.dat / .out)
        if strcmp(data_format, 'vd')
            noise_raw = load_noise_file_vd(noise_filepath, n_coils);
        else
            noise_raw = load_noise_file_legacy(noise_filepath, n_coils);
        end
        fprintf('  Source : separate file  (%d samples x %d coils)\n', ...
            size(noise_raw,1), size(noise_raw,2));
end

% Fallback (or explicit selection): outer k-space lines from signal scan
if strcmp(noise_source, 'kspace')
    % Use outer fraction of ky lines where signal power is negligible.
    % kspace_sig is [Nx, Ny_full, Ncoil] from the load function.
    Nky    = size(kspace_sig, 2);
    n_edge = max(1, round(kspace_noise_fraction * Nky));
    idx_lo = 1:n_edge;
    idx_hi = (Nky - n_edge + 1):Nky;
    noise_k = kspace_sig(:, [idx_lo, idx_hi], :);   % [Nx, 2*n_edge, Ncoil]
    noise_raw = reshape(noise_k, [], n_coils);        % [Nsamples x Ncoil]
    fprintf('  Source : outer k-space  (%d lines, %d samples x %d coils)\n', ...
        2*n_edge, size(noise_raw,1), size(noise_raw,2));
end

% Strip dummy channel from noise now that noise_raw is defined.
% The noise file also has 9 channels; remove channel 1 to match images.
if dummy_ch && size(noise_raw, 2) == 9
    noise_raw = noise_raw(:, 2:end);
    fprintf('  9-channel noise: channel 1 (dummy) removed\n');
end

if size(noise_raw, 2) ~= n_coils
    error('Noise data has %d coils but signal has %d.', size(noise_raw,2), n_coils);
end

%% ========================
%  NOISE COVARIANCE
%% ========================
fprintf('Computing noise covariance...\n');
Psi = (noise_raw' * noise_raw) / size(noise_raw, 1);   % [Ncoil x Ncoil]

if prewhiten
    L_full = chol(Psi, 'lower');
    fprintf('  Pre-whitening enabled.  Cholesky cond: %.1f\n\n', cond(L_full));
else
    L_full = eye(n_coils);
    fprintf('  Pre-whitening disabled.\n\n');
end

% Full whitened image matrix (used for coil selection if n_coils=8)
A_full   = reshape(images, [], n_coils);    % [Nvox x Ncoil]
A_full_w = A_full / L_full';                % whitened

%% ========================
%  COIL SUBSET SELECTION  (8-channel only)
%% ========================
if n_coils == 8
    if ~isempty(manual_channels)
        % ?? Manual override ??????????????????????????????????????????????
        if numel(manual_channels) ~= 4 || any(manual_channels < 1) || any(manual_channels > 8)
            error('manual_channels must be a 1x4 vector with values in [1,8]. Got: %s', ...
                mat2str(manual_channels));
        end
        selected_ch = manual_channels(:)';
        fprintf('Channel selection: MANUAL  channels=%s\n\n', mat2str(selected_ch));

    else
        % ?? Automatic selection ??????????????????????????????????????????
        fprintf('Selecting optimal 4-channel subset...\n');

    % Candidate subsets: alternating (encircling) and neighbouring
    subsets = {
        [1,3,5,7], [2,4,6,8], ...               % alternating
        [1,2,3,4],[2,3,4,5],[3,4,5,6], ...       % neighbouring
        [4,5,6,7],[5,6,7,8],[6,7,8,1], ...
        [7,8,1,2],[8,1,2,3]                       % wrap-around
    };
    subset_labels = {'Alt-1357','Alt-2468', ...
        'Nbr-1234','Nbr-2345','Nbr-3456', ...
        'Nbr-4567','Nbr-5678','Nbr-6781', ...
        'Nbr-7812','Nbr-8123'};

    % Two range values:
    %   hw_range_full : actual hardware span ? used for subset SELECTION.
    %                   No buffer here; the buffer is a hardware-setting
    %                   safety margin, not a measure of achievable range.
    %   hw_range_buf  : buffered span ? used only for final feasibility
    %                   warnings in the phase rotation block below.
    hw_range_full = hw_ph_hi - hw_ph_lo;
    hw_range_buf  = hw_range_full - 2*hw_ph_buf;

    % SNR mask: use full-whitened data (noise-normalised, all 8 channels)
    snr_quick  = sqrt(sum(abs(A_full_w).^2, 2));
    mask_quick = snr_quick > snr_threshold * max(snr_quick);

    % Raw (un-whitened) full data for per-subset whitening in the loop.
    % Each subset must be whitened by its OWN 4-channel noise sub-block,
    % not by the full 8-channel Cholesky factor.  Using A_full_w columns
    % with a 4-channel un-whitening step is inconsistent: the full
    % whitening absorbs cross-channel correlations from all 8 channels,
    % including those outside the subset, producing wrong phase estimates
    % and artificially large score differences between subsets.
    A_full_raw = reshape(images, [], n_coils);   % [Nvox x 8] un-whitened

    scores  = zeros(1, numel(subsets));
    margins = zeros(numel(subsets), 4);

    fprintf('  %-12s  Score   Mode-1   Mode-2   Mode-3   Mode-4\n', 'Subset');
    fprintf('  %s\n', repmat('-', 1, 58));

    for s = 1:numel(subsets)
        ch = subsets{s};

        % Per-subset whitening: extract raw channels then whiten with
        % the 4x4 noise sub-block for this channel combination only.
        A_raw_s = A_full_raw(mask_quick, ch);    % [Nvox_mask x 4] raw
        L_s     = chol(Psi(ch,ch), 'lower');     % 4x4 subset Cholesky
        A_s_w   = A_raw_s / L_s';               % properly whitened subset

        [~, ~, V_s] = svd(A_s_w, 'econ');
        C_s = L_s' \ V_s;                      % [4 channels x 4 modes], C_s(:,m) = mode m

        for m = 1:4
            phases_m = rad2deg(angle(C_s(:,m)));   % column m = mode m weights across channels
            [~, arc_m, ~] = min_circular_arc(phases_m);
            margins(s,m) = hw_range_full - arc_m;
            scores(s) = scores(s) + mode_weights(m) * max(0, margins(s,m));
        end

        fprintf('  %-12s  %5.1f  %+7.1f  %+7.1f  %+7.1f  %+7.1f\n', ...
            subset_labels{s}, scores(s), ...
            margins(s,1), margins(s,2), margins(s,3), margins(s,4));
    end

    % Three-tier selection (Mode 1 weighted above Mode 2):
    %   Tier 3 : both Mode 1 and Mode 2 fit in full hardware range  (best)
    %   Tier 2 : Mode 1 fits, Mode 2 does not
    %   Tier 1 : Mode 2 fits, Mode 1 does not
    %   Tier 0 : neither fits
    % Margin-weighted score breaks ties within each tier.
    tier = (margins(:,1) > 0) * 2 + (margins(:,2) > 0);
    composite = tier * 1e6 + scores(:);
    [~, best] = max(composite);

    selected_ch = subsets{best};
    tier_labels = {'neither feasible','Mode 2 only','Mode 1 only','Modes 1 & 2'};
    fprintf('\n  => Best subset: %s  channels=%s  (%s in hardware range)\n\n', ...
        subset_labels{best}, mat2str(selected_ch), tier_labels{tier(best)+1});

    end   % end: automatic selection

else
    selected_ch = 1:4;
    fprintf('4-channel data ? using all channels.\n\n');
end

% Extract selected channels and their noise covariance
images_sel = images(:, :, :, selected_ch);
Psi_sel    = Psi(selected_ch, selected_ch);

if prewhiten
    L_sel = chol(Psi_sel, 'lower');
else
    L_sel = eye(4);
end

A_sel   = reshape(images_sel, [], 4);
A_sel_w = A_sel / L_sel';

%% ========================
%  ROI MASK
%% ========================
fprintf('Building ROI mask (type: %s)...\n', roi_type);

snr_w_vol = reshape(sqrt(sum(abs(A_sel_w).^2, 2)), Nx, Ny, Nz);

% SNR-based mask
mask_snr = snr_w_vol > snr_threshold * max(snr_w_vol(:));

% Radial mask: keep inner radial_fraction * R per slice
mask_rad = true(Nx, Ny, Nz);
if strcmp(roi_type,'radial') || strcmp(roi_type,'both')
    [cx, cy] = phantom_center(snr_w_vol(:,:,ceil(Nz/2)));
    [XX, YY] = meshgrid(1:Ny, 1:Nx);
    R_est    = min([cx-1, Nx-cx, cy-1, Ny-cy]);
    R_map    = sqrt((XX-cy).^2 + (YY-cx).^2);
    mask_2d  = R_map <= radial_fraction * R_est;
    mask_rad = repmat(mask_2d, 1, 1, Nz);
end

mask_manual = false(Nx,Ny,Nz);

if strcmp(roi_type,'manual') && ~load_roi

    sl = ceil(Nz/2);

    figure;
    imagesc(abs(snr_w_vol(:,:,sl)));
    axis image;
    colormap gray;
    colorbar;
    title('Draw ROI and double-click inside when finished');

    h = drawcircle('Color','r');
    wait(h);

    roi2d = createMask(h);

    close(gcf)

    % Apply same ROI to all slices
    mask_manual = repmat(roi2d,[1 1 Nz]);

    B = bwboundaries(mask_manual);

    % % helper function
    % plot_roi = @(B) cellfun(@(b) ...
    %     plot(b(:,2), b(:,1), 'r', 'LineWidth', 1.5), ...
    %     B, 'UniformOutput', false);

   

end

if load_roi

    fprintf('Loading ROI from %s\n', roi_filename);

    S = load(roi_filename);

    if ~isfield(S,'mask')
        error('ROI file does not contain variable ''mask''.');
    end

    mask = logical(S.mask);

     B = bwboundaries(mask);

    % % helper function
    % plot_roi = @(B) cellfun(@(b) ...
    %     plot(b(:,2), b(:,1), 'r', 'LineWidth', 1.5), ...
    %     B, 'UniformOutput', false);

    if Nz == 1
        % 2D dataset
        if ~isequal(size(mask),[Nx Ny])
            error('ROI dimensions do not match current dataset.');
        end
    elseif Nz>1
        if ~isequal(size(mask),[Nx Ny Nz])
            error('ROI dimensions do not match current dataset.');
        end
    end
    fprintf('  Loaded ROI: %d voxels\n',nnz(mask));

elseif ~load_roi

    switch roi_type
        case 'snr'
            mask = mask_snr;

        case 'radial'
            mask = mask_rad;

        case 'both'
            mask = mask_snr & mask_rad;

        case 'manual'
            mask = mask_manual;

            if save_roi

                save(roi_filename,'mask');

                fprintf('ROI saved to %s\n',roi_filename);

            end

        otherwise
            error('Unknown roi_type: %s', roi_type);
    end
end

fprintf('  %d / %d voxels in ROI (%.1f%%)\n\n', ...
    sum(mask(:)), numel(mask), 100*sum(mask(:))/numel(mask));

%% ========================
%  PCA / SVD
%% ========================
fprintf('Computing compression matrix via SVD...\n');
A_roi = A_sel_w(mask(:), :);   % [Nvox_roi x 4], whitened + masked

[~, S_svd, V_w] = svd(A_roi, 'econ');   % V_w: [4 x 4]

sv      = diag(S_svd);
sv_norm = sv / sv(1);
snr2_cum = cumsum(sv.^2) / sum(sv.^2);

fprintf('  Mode  Norm-SV  SNR²-contrib  Cumulative\n');
for m = 1:4
    fprintf('    %d   %.4f     %5.1f%%        %5.1f%%\n', ...
        m, sv_norm(m), 100*(sv(m)^2)/sum(sv.^2), 100*snr2_cum(m));
end
fprintf('\n');

% Compression matrix in original (non-whitened) space.
% Convention matches calculate_settings_complex: A(channel, mode).
% V_w: [4 channels x 4 modes], modes in columns.
% C_orig = V_w / L_sel' maps whitened singular vectors back to the
% original channel space, preserving [channels x modes] orientation.
C_orig = L_sel' \ V_w;    % [4 channels x 4 modes], C_orig(:,m) = mode m weights

%% ========================
%  PHASE FEASIBILITY AND PER-MODE ROTATION
%
% For each mode row: find the optimal global phase rotation that maps all
% 4 element phases into the hardware buffered range [ph_lo_buf, ph_hi_buf].
%
% Key insight: a global rotation of mode row m by e^{j*theta} does not
% change the spatial sensitivity pattern of the mode ? it only shifts its
% phase reference.  This degree of freedom is used to fit the mode into
% the hardware range without sacrificing imaging performance.
%
% A mode is infeasible if its minimum circular arc > hardware range.
% This is independent of rotation; rotation cannot reduce the span.
%% ========================
fprintf('Phase feasibility and rotation:\n');
fprintf('  Hardware buffered range: [%+.1f, %+.1f] deg  (span %.1f deg)\n\n', ...
    hw_ph_lo + hw_ph_buf, hw_ph_hi - hw_ph_buf, ...
    (hw_ph_hi - hw_ph_lo) - 2*hw_ph_buf);

ph_lo_buf  = hw_ph_lo + hw_ph_buf;
ph_hi_buf  = hw_ph_hi - hw_ph_buf;
ph_ctr     = (ph_lo_buf + ph_hi_buf) / 2;
hw_rng_buf = ph_hi_buf - ph_lo_buf;

A_out = complex(zeros(4));  % [4 channels x 4 modes]: A_out(channel,mode)
                             % column m = mode m weights across physical channels

for m = 1:4
    col_m  = C_orig(:,m);                   % column m = mode m weights [4 channels x 1]
    ph_m   = rad2deg(angle(col_m));

    [~, arc_deg, arc_start_deg] = min_circular_arc(ph_m);
    arc_ctr_deg = arc_start_deg + arc_deg/2;
    theta_deg   = ph_ctr - arc_ctr_deg;

    % Three feasibility levels:
    %   'ok'      : fits within buffered range ? full safety margin
    %   'tight'   : fits in full hardware range but within buffer zone
    %   'bad'     : exceeds full hardware range ? cannot be implemented
    if arc_deg <= hw_rng_buf
        level      = 'ok';
        status_str = sprintf('OK      span %5.1f deg  buffer margin %+.1f deg', ...
            arc_deg, hw_rng_buf - arc_deg);
    elseif arc_deg <= (hw_ph_hi - hw_ph_lo)
        level      = 'tight';
        status_str = sprintf('TIGHT   span %5.1f deg  in hardware range, within buffer zone', ...
            arc_deg);
        if m <= 2
            fprintf('    Note: Mode %d is within full hardware range but inside the %.0f deg buffer zone.\n', ...
                m, hw_ph_buf);
            fprintf('    Hardware can implement it; boundary calibration accuracy may be reduced.\n');
        end
    else
        level      = 'bad';
        status_str = sprintf('INFEASIBLE  span %.1f deg exceeds hardware range %.1f deg', ...
            arc_deg, hw_ph_hi - hw_ph_lo);
        if m <= 2
            warning('Mode %d exceeds the full hardware phase range and cannot be correctly implemented.', m);
        else
            fprintf('    (Mode %d infeasibility has lower impact ? consider 2-mode compression.)\n', m);
        end
    end

    col_rot    = col_m * exp(1j * deg2rad(theta_deg));
    A_out(:,m) = col_rot;                   % assign rotated weights to column m

    ph_rot = rad2deg(angle(col_rot));
    fprintf('  Mode %d  [%+6.1f %+6.1f %+6.1f %+6.1f] deg  %s\n', ...
        m, ph_rot(1), ph_rot(2), ph_rot(3), ph_rot(4), status_str);
end
fprintf('\n');

%% ========================
%  SNR RETENTION ANALYSIS
%
% Reference: optimal matched-filter combination of all 4 whitened channels
%   SNR_ref(x)  = ||a_w(x)||_2
%
% Compressed M modes:
%   SNR_M(x) = ||V_w(:,1:M)' * a_w(x)||_2
%   Retention = SNR_M(x) / SNR_ref(x)
%
% Note: V_w are the right singular vectors of the whitened data matrix.
% Since noise is identity in the whitened domain, this is the correct
% noise-normalised SNR rather than a signal-only norm.
%% ========================
fprintf('SNR retention analysis:\n');

snr2_ref = sum(abs(A_sel_w).^2, 2) + eps;   % [Nvox x 1]

combos      = {1:4,  1:3,  1:2,  1};
combo_names = {'All 4 modes','Modes 1-3','Modes 1-2','Mode 1 only'};
ret_maps    = cell(1,4);

fprintf('  %-14s  Mean    Median   5th-pct   Min\n', 'Combination');
fprintf('  %s\n', repmat('-',1,54));

for k = 1:4
    M       = combos{k};
    proj    = A_sel_w * V_w(:,M);            % [Nvox x |M|]
    snr2_M  = sum(abs(proj).^2, 2);
    ret     = sqrt(snr2_M ./ snr2_ref);      % [Nvox x 1], SNR ratio
    ret_maps{k} = reshape(ret, Nx, Ny, Nz);

    r_roi = ret(mask(:));   % statistics within ROI only
    fprintf('  %-14s  %.3f   %.3f    %.3f     %.3f\n', ...
        combo_names{k}, mean(r_roi), median(r_roi), ...
        prctile(r_roi,5), min(r_roi));
end
fprintf('\n');

%% ========================
%  VISUALISATION
%% ========================
if plot_results
    sl = ceil(Nz/2);   % central slice for 2-D display

    % --- Singular value spectrum ---
    figure('Name','SVD Spectrum','Position',[100 100 700 300]);
    subplot(1,2,1);
    bar(sv_norm, 'FaceColor',[0.2 0.4 0.8]);
    set(gca,'XTick',1:4,'XTickLabel',{'Mode 1','Mode 2','Mode 3','Mode 4'});
    ylabel('Normalised singular value'); title('SVD spectrum');

    subplot(1,2,2);        hold on
    bar(100*(sv.^2)/sum(sv.^2), 'FaceColor',[0.8 0.3 0.2]);
    set(gca,'XTick',1:4,'XTickLabel',{'Mode 1','Mode 2','Mode 3','Mode 4'});
    ylabel('SNR² contribution (%)'); title('Per-mode SNR² contribution');

    % --- SNR retention maps ---
    figure('Name','SNR Retention','Position',[100 450 1100 280]);
    for k = 1:4
        subplot(1,4,k);
        imagesc(ret_maps{k}(:,:,sl), [0 1]);
        colormap(hot); colorbar; axis image off;
        title(combo_names{k},'FontSize',9);
        
        if strcmp(roi_type,'manual')
            hold on
            plot_roi(B)
            hold off
        end
    end
    sgtitle('SNR retention vs optimal 4-channel matched filter');

    % --- Virtual mode sensitivity maps (magnitude) ---
    figure('Name','Mode Maps','Position',[100 750 1100 280]);
    for m = 1:4
        mode_img = reshape(A_sel_w * V_w(:,m), Nx, Ny, Nz);
        subplot(1,4,m);
        imagesc(abs(mode_img(:,:,sl)));
        colormap(gray); colorbar; axis image off;
        title(sprintf('Mode %d  (SV=%.3f)', m, sv_norm(m)), 'FontSize',9);
        if strcmp(roi_type,'manual')
            hold on
            plot_roi(B)
            hold off
        end
    end
    sgtitle('Virtual mode sensitivity maps ? magnitude (whitened domain)');

    % --- Diagnostic: complex-domain mode maps ---
    % If a mode appears as pure noise in the magnitude display but has
    % spatial structure in its real or imaginary part, the cause is complex
    % cancellation (the mode is valid but the channels partially cancel in
    % both Re and Im, leaving near-zero magnitude everywhere).
    % If Re, Im, and magnitude all look like noise, the mode itself is
    % dominated by noise ? likely a dead/anomalous channel.
    figure('Name','Mode Maps ? Complex Diagnostic','Position',[100 1060 1600 520]);
    for m = 1:4
        mode_img = reshape(A_sel_w * V_w(:,m), Nx, Ny, Nz);
        sl_img   = mode_img(:,:,sl);

        subplot(3,4,m);
        imagesc(abs(sl_img));   colormap(gray); colorbar; axis image off;
        title(sprintf('Mode %d  |z|  SV=%.3f', m, sv_norm(m)), 'FontSize',8);
       
        if strcmp(roi_type,'manual')
            hold on
            plot_roi(B)
            hold off
        end

        subplot(3,4,m+4);
        imagesc(real(sl_img));  colormap(gray); colorbar; axis image off;
        title(sprintf('Mode %d  Re(z)', m), 'FontSize',8);
        if strcmp(roi_type,'manual')
            hold on
            plot_roi(B)
            hold off
        end
        subplot(3,4,m+8);
        imagesc(imag(sl_img));  colormap(gray); colorbar; axis image off;
        title(sprintf('Mode %d  Im(z)', m), 'FontSize',8);
        if strcmp(roi_type,'manual')
            hold on
            plot_roi(B)
            hold off
        end
    end
    sgtitle('Complex diagnostic: |z|, Re(z), Im(z) per mode');

    % --- Diagnostic: individual whitened channel images ---
    % Reveals dead channels, phase anomalies, or data loading errors
    % before they propagate into the SVD.
    figure('Name','Channel Images (whitened)','Position',[100 100 1100 540]);
    for ch = 1:4
        % Magnitude
        ch_img = reshape(A_sel_w(:,ch), Nx, Ny, Nz);
        subplot(2,4,ch);
        imagesc(abs(ch_img(:,:,sl)));
        colormap(gray); colorbar; axis image off;
        title(sprintf('Ch %d (orig %d) |z|', ch, selected_ch(ch)), 'FontSize',8);
        if strcmp(roi_type,'manual')
            hold on
            plot_roi(B)
            hold off
        end
        % Phase
        subplot(2,4,ch+4);
        imagesc(angle(ch_img(:,:,sl)), [-pi pi]);
        colormap(hsv); colorbar; axis image off;
        title(sprintf('Ch %d  phase', ch), 'FontSize',8);
        if strcmp(roi_type,'manual')
            hold on
            plot_roi(B)
            hold off
        end
    end
    sgtitle('Individual whitened channel images ? magnitude and phase');

    % --- Noise correlation matrix ---
    figure('Name','Noise Correlation','Position',[820 100 350 300]);
    Psi_norm = abs(Psi_sel) ./ sqrt(diag(Psi_sel)*diag(Psi_sel)');
    imagesc(Psi_norm, [0 1]); colormap(hot); colorbar; axis square;
    title('Noise correlation matrix'); xlabel('Channel'); ylabel('Channel');
    set(gca,'XTick',1:4,'YTick',1:4);
end

%% ========================
%  SAVE OUTPUT
%% ========================
if save_compression_matrix
    A = A_out;   % [4 x 4] complex, phase-rotated to fit hardware range

    % Build filename: base name + channel tag, saved in signal file folder.
    % e.g. compression_matrix_ch3456.mat
    ch_tag      = sprintf('ch%s', sprintf('%d', selected_ch));
    [~, base, ~] = fileparts(output_file);
    output_path  = fullfile(sig_path, sprintf('%s_%s.mat', base, ch_tag));

    save(output_path, 'A', ...
        'sv', 'sv_norm', 'snr2_cum', ...       % SVD diagnostics
        'selected_ch', ...                      % which physical channels used
        'prewhiten', 'roi_type', ...            % pipeline settings
        'ret_maps', 'combo_names');             % SNR retention maps

    fprintf('Compression matrix saved: %s\n', output_path);
    fprintf('  Load with loadC=1 in calculate_settings_complex.m\n\n');
end
%% ========================
%  LOCAL FUNCTIONS
%% ========================

% ?? RAW DATA LOADERS ?????????????????????????????????????????????????????

function [images, kspace_out, params] = load_pca_signal_vd(filepath, params, config)
% Load VD13+ format signal scan via mapVBVD.
% Returns multi-channel images [Nx, Ny, Nz, Ncoil] and raw k-space.
%
% Unlike the B1+ pipeline, which reads multiple files at different flip
% angles, the PCA pipeline reads a single signal file and preserves all
% channels separately.  The raw k-space is also returned so the caller
% can extract noise from its outer lines if needed.

    twix = mapVBVD(filepath);

    twix.image.flagRemoveOS = config.remove_oversampling;
    kspace = single(twix.image());   % dimension order: Col, Cha, Lin, Par, Sli, ...

    fprintf('    Raw k-space: %s\n', mat2str(size(kspace)));

    % Average over acquisition averages (mapVBVD dim 6 = Ave) before squeeze.
    % size(x, d) returns 1 safely when x has fewer than d dimensions, so
    % no ndims guard is needed.
    if size(kspace, 6) > 1
        n_avg_ks = size(kspace, 6);
        kspace   = mean(kspace, 6);
        fprintf('    Averaged %d acquisitions (k-space)\n', n_avg_ks);
    end
    kspace = squeeze(kspace);

    % Standardise to [Nx, Ny, Ncoil] or [Nx, Ny, Nz, Ncoil]
    % mapVBVD typical layout after squeeze: Col × Cha × Lin [× Par × Sli]
    ndim = ndims(kspace);
    sz   = size(kspace);

    switch ndim
        case 3
            % [Col, Cha, Lin] ? permute to [Col, Lin, Cha] = [Nx, Ny, Ncoil]
            kspace = permute(kspace, [1, 3, 2]);
            fprintf('    Layout: [Nx=%d, Ny=%d, Ncoil=%d]\n', sz(1), sz(3), sz(2));

        case 4
            % [Col, Cha, Lin, Sli] ? [Nx, Ny, Ncoil, Nz] ? [Nx, Ny, Nz, Ncoil]
            kspace = permute(kspace, [1, 3, 4, 2]);
            fprintf('    Layout: [Nx=%d, Ny=%d, Nz=%d, Ncoil=%d]\n', ...
                sz(1), sz(3), sz(4), sz(2));

        otherwise
            error('Unexpected k-space dimensionality: %d', ndim);
    end

    % Save raw k-space before FFT (used for outer-line noise extraction)
    kspace_out = double(kspace);

    % 2-D FFT along first two dimensions (readout x phase-encode)
    images = fftshift(fftshift( ...
        ifft(ifft(kspace, [], 1), [], 2), 1), 2);
    images = double(images);

    % Extract scan parameters
    try
        header = eval_twix_hdr(filepath);
        params.flip_angle   = get_flip_angle(header);
        params.frequency    = get_field_safe(header, ...
            {'Dicom.lFrequency','Meas.lFrequency'}, 127.7e6);
        params.field_strength = get_field_safe(header, ...
            {'Dicom.flMagneticFieldStrength','Meas.flMagneticFieldStrength'}, 3.0);
    catch
        % Non-critical; continue without header parameters
    end
end

function [images, kspace_out, params] = load_pca_signal_legacy(filepath, params)
% Load legacy format signal scan via read_meas_dat_UNIX.

    opt.ReturnStruct = 1;
    meas   = read_meas_dat_UNIX(filepath, opt);
    images = double(mrir_conventional_2d(meas.data));

    % Legacy reader does not easily expose raw k-space; return empty.
    % Outer-k-space noise extraction is not available for legacy format.
    kspace_out = [];

    try
        params.flip_angle     = meas.prot.adFlipAngleDegree;
        params.frequency      = meas.prot.lFrequency;
        params.field_strength = meas.prot.flNominalB0;
    catch
    end
end

function noise_raw = load_noise_file_vd(filepath, n_coils_expected)
% Load a separate noise-only .dat file in VD format.
% The file typically contains a short noise-only scan (RF off).

    twix = mapVBVD(filepath);

    % Try twix.noise() first, then twix.image() as fallback
    if isfield(twix, 'noise') && ~isempty(twix.noise)
        raw = double(squeeze(twix.noise()));
    else
        raw = double(squeeze(twix.image()));
    end

    % Ensure [Nsamples x Ncoil]
    if ndims(raw) == 3
        raw = reshape(raw, [], size(raw,2));   % flatten readout × lines
    end
    if size(raw,2) ~= n_coils_expected && size(raw,1) == n_coils_expected
        raw = raw.';
    end

    noise_raw = raw;
end

function noise_raw = load_noise_file_legacy(filepath, n_coils_expected)
% Load a separate noise-only .dat file in legacy format.

    opt.ReturnStruct = 1;
    meas     = read_meas_dat_UNIX(filepath, opt);
    noise_raw = reshape(double(meas.data), [], n_coils_expected);
end

function [feasible, arc_deg, arc_start_deg] = min_circular_arc(phases_deg)
% Minimum circular arc containing all phase values.
%
% Uses the largest-gap algorithm:
%   sort phases on [0,360); the minimum enclosing arc is
%   360 deg minus the largest gap between consecutive phases.
%
% phases_deg   : 1 x N  phase values in degrees (any range)
% feasible     : true if min arc <= hw_range_deg (set by caller)
% arc_deg      : width of the minimum enclosing arc in degrees
% arc_start_deg: angle at which the minimum arc begins

    ph = mod(phases_deg(:)', 360);   % normalise to [0, 360)
    ph_s = sort(ph);
    n    = numel(ph_s);

    gaps = [diff(ph_s), 360 - ph_s(end) + ph_s(1)];
    [max_gap, idx] = max(gaps);

    arc_deg       = 360 - max_gap;
    arc_start_deg = ph_s(mod(idx, n) + 1);   % phase after the largest gap

    % feasibility is evaluated by the caller against hw_range_buf
    feasible = true;   % placeholder; caller compares arc_deg to threshold
end

function [cx, cy] = phantom_center(snr_slice)
% Intensity-weighted centroid of a 2-D SNR map.
    thresh   = snr_slice > 0.2 * max(snr_slice(:));
    [rr, cc] = find(thresh);
    w        = snr_slice(thresh);
    w        = w / sum(w);
    cx = round(dot(rr, w));
    cy = round(dot(cc, w));
end

function plot_roi(B)

    for k = 1:numel(B)
        b = B{k};
        plot(b(:,2), b(:,1), ...
            'r', 'LineWidth', 1.5);
    end

end