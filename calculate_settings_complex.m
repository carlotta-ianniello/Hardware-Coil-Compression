% =========================================================
% calculate_settings_complex.m
%
% Programmable mixing matrix ? target-to-hardware solver.
% Operates entirely in linear complex domain.
%
% Key design principles:
%   - Cal_mag + Cal_ph converted to Cal_complex immediately after load;
%     no dB or phase arithmetic thereafter until output display.
%   - Solver: Euclidean distance in complex plane ? no weights to tune,
%     phase error automatically down-weighted at high attenuation.
%   - Amplitude scaling: per-element-aware.  The global scale is set by
%     the single most-constrained active element relative to its own
%     per-element calibration ceiling ? not by global_max_amp applied
%     uniformly.  This uses the full available calibration range.
%   - amp_buffer defaults to 0.  Set to a small positive value (e.g.
%     0.02, ?0.18 dB) only if you need a margin against calibration drift.
%   - Phase boundary: symmetric buffer at both ends of the intersection
%     of achievable phase ranges across all N×N elements.
%   - dB and degrees appear only at input and output display.
% =========================================================
clear

%% =========================
% EXPORT FLAG
%% =========================
export_to_IDE = 0;
dothfilename  = 'C:\Users\iannic01\Personal\Research\Arduino_test\modular_board_test\4x4_matrix\4x4_matrix_settings\PCA7812_amp_ROI3';

%% =========================
% MODE FLAGS
%% =========================
loadC     = 1;
randC     = 0;
useCalVal = 0;

if loadC + randC + useCalVal > 1
    error('Select at most one of: loadC, randC, useCalVal');
end

%% =========================
% HARDWARE PARAMETERS
%% =========================
N         = 4;
n_att     = 64;
n_ph      = 13;
n_interp  = 64;
att_step  = 0.5;      % dB per attenuator LSB
V_max     = 12;       % maximum phase shifter voltage (V)

% Amplitude buffer (fractional, linear).
% Default 0: targets reach the exact per-element calibration boundary.
% Set e.g. 0.02 (~0.18 dB) to guard against calibration drift.
amp_buffer = 0;

% Symmetric phase buffer, equal margin at each end of the common range.
ph_buffer_deg = 10;

% Phase reporting threshold.
ph_report_thresh_dB = -30;

%% =========================
% PATH + LOAD CALIBRATION
%% =========================
cal_path = 'C:\Users\iannic01\OneDrive - NYU Langone Health\4x4 matrix settings calculation/';
addpath(genpath(cal_path));
load([cal_path 'calibration_4x4.mat']);

%% =========================
% DIAGNOSTIC: Downsample 64×13 calibration entries to 64×7
% by removing even-indexed columns (1 V, 3 V, 5 V, 7 V, 9 V, 11 V steps),
% retaining columns 1,3,5,7,9,11,13 ? 0, 2, 4, 6, 8, 10, 12 V.
%
% Run AFTER loading calibration, BEFORE the conversion to Cal_complex.
% Set downsample_test = 1 to enable; 0 to use the original calibration.
%% =========================
downsample_test = 0;

if downsample_test
    cols_keep = 1:2:13;   % [1, 3, 5, 7, 9, 11, 13]
    n_downsampled = 0;
    for row = 1:N
        for col = 1:N
            if size(calibration(row,col).C_mag, 2) == 13
                calibration(row,col).C_mag = calibration(row,col).C_mag(:, cols_keep);
                calibration(row,col).C_ph  = calibration(row,col).C_ph(:,  cols_keep);
                n_downsampled = n_downsampled + 1;
            end
        end
    end
    fprintf('Downsampled %d element(s) from 64×13 to 64×7.\n\n', n_downsampled);
end
%% =========================
% CONVERT CALIBRATION TO LINEAR COMPLEX
%
% Interpolation axes are derived from each element's ACTUAL dimensions
% rather than the global n_att / n_ph parameters.  This makes the code
% robust to calibration matrices with inconsistent sizes across elements
% (e.g. an incomplete calibration run that produced fewer voltage steps
% for some elements).
%
% A diagnostic table is printed so any size mismatch is immediately
% visible.  All elements are still interpolated to the common n_interp
% grid so downstream code is unaffected.
%% =========================
Cal_complex = cell(N, N);

% --- Diagnostic: print actual calibration sizes ---
fprintf('Calibration matrix dimensions [n_att x n_ph]:\n');
expected_sz = [];
size_mismatch = false;
for row = 1:N
    fprintf('  ');
    for col = 1:N
        sz = size(calibration(row,col).C_mag);
        fprintf('(%d,%d):%dx%d  ', row, col, sz(1), sz(2));
        if isempty(expected_sz)
            expected_sz = sz;
        elseif ~isequal(sz, expected_sz)
            size_mismatch = true;
        end
    end
    fprintf('\n');
end
if size_mismatch
    warning(['Calibration elements have inconsistent sizes (shown above). ' ...
             'Each element will be interpolated from its own grid. ' ...
             'Consider recalibrating the affected elements.']);
else
    fprintf('  All elements: %d x %d  (consistent)\n', expected_sz(1), expected_sz(2));
end
fprintf('\n');

% --- Conversion loop ---
for row = 1:N
    for col = 1:N

        % Actual dimensions for this element
        [n_att_el, n_ph_el] = size(calibration(row,col).C_mag);

        % Per-element interpolation axes derived from actual size
        v_old_el = linspace(1, n_ph_el,  n_ph_el);
        v_new_el = linspace(1, n_ph_el,  n_interp);
        a_old_el = linspace(1, n_att_el, n_att_el);
        a_new_el = linspace(1, n_att_el, n_interp);

        % Single-line dB/degrees -> linear complex
        raw = 10.^(calibration(row,col).C_mag / 20) ...
           .* exp(1j * deg2rad(calibration(row,col).C_ph));  % [n_att_el x n_ph_el]

        % Interpolate voltage axis (columns: n_ph_el -> n_interp)
        raw_i = interp1(v_old_el, real(raw).', v_new_el, 'pchip').' ...
             + 1j*interp1(v_old_el, imag(raw).', v_new_el, 'pchip').'; % [n_att_el x n_interp]

        % Interpolate attenuator axis (rows: n_att_el -> n_interp)
        Cal_complex{row,col} = ...
              interp1(a_old_el, real(raw_i), a_new_el, 'pchip') ...
           + 1j*interp1(a_old_el, imag(raw_i), a_new_el, 'pchip');    % [n_interp x n_interp]
    end
end

%% =========================
% AMPLITUDE BOUNDARIES  (per-element, from ORIGINAL calibration data)
%
% max_amp and min_amp are computed from the raw calibration measurements,
% NOT from the interpolated Cal_complex surface.
%
% Rationale: pchip interpolation of Re and Im components separately can
% produce magnitudes abs(Re + j*Im) that exceed any value present in the
% original data (because each component is bounded but the vector magnitude
% is not).  Using the interpolated surface can therefore report a max_amp
% that is physically unachievable, causing targets to be set above the
% true hardware ceiling.  Using the original data guarantees every boundary
% value corresponds to a real calibration measurement.
%
% Additionally: max_amp must occur at i=1 (minimum attenuation).  A warning
% is issued if this is violated, which would indicate a calibration error.
%% =========================
max_amp = zeros(N, N);
min_amp = zeros(N, N);

for row = 1:N
    for col = 1:N
        % Original (pre-interpolation) amplitudes
        raw_amp = 10.^(calibration(row,col).C_mag / 20);   % [n_att x n_ph]

        max_amp(row,col) = max(raw_amp(:));
        min_amp(row,col) = min(raw_amp(:));

        % Sanity check: max should occur at row 1 (minimum attenuation).
        % If it occurs elsewhere, the calibration data may have errors.
        [~, max_idx] = max(raw_amp(:));
        [i_at_max, ~] = ind2sub(size(raw_amp), max_idx);
        if i_at_max ~= 1
            warning(['Element (%d,%d): max amplitude occurs at attenuator ' ...
                     'index %d (%.1f dB), not index 1 (0 dB). ' ...
                     'Calibration data may contain errors.'], ...
                row, col, i_at_max, (i_at_max-1)*att_step);
        end
    end
end

global_max_amp = min(max_amp(:));
global_min_amp = max(min_amp(:));

fprintf('Per-element amplitude ceilings (dB):\n');
for row = 1:N
    for col = 1:N
        fprintf('  (%d,%d)  max %.2f dB   min %.2f dB\n', row, col, ...
            20*log10(max_amp(row,col)), 20*log10(min_amp(row,col)));
    end
end
fprintf('Global limits:  max %.2f dB   min %.2f dB\n\n', ...
    20*log10(global_max_amp), 20*log10(global_min_amp));

%% =========================
% PHASE BOUNDARIES  (symmetric intersection)
%% =========================
ph_min_rad = zeros(N, N);
ph_max_rad = zeros(N, N);

for row = 1:N
    for col = 1:N
        % Use ORIGINAL calibration phase at row 1 (minimum attenuation),
        % not the interpolated surface, for the same reason as max_amp:
        % pchip interpolation of Re/Im can create phase extrema outside
        % the true hardware range.
        phases_orig = angle( 10.^(calibration(row,col).C_mag(1,:)/20) ...
                          .* exp(1j*deg2rad(calibration(row,col).C_ph(1,:))) );  % 1 x n_ph_el

        % Detect actual wrap-around via discontinuity in original data.
        % Large total span (> 180 deg) is expected for JSPHS-42+ (~204 deg)
        % and is NOT an error.
        if any(abs(diff(phases_orig)) > deg2rad(170))
            warning(['Element (%d,%d): possible phase wrap-around in original ' ...
                     'calibration data (jump > 170 deg between adjacent voltage ' ...
                     'points).  Phase boundary may be incorrect.'], row, col);
        end

        ph_min_rad(row,col) = min(phases_orig);
        ph_max_rad(row,col) = max(phases_orig);
    end
end

global_ph_lo  = max(ph_min_rad(:));
global_ph_hi  = min(ph_max_rad(:));
ph_range_rad  = global_ph_hi - global_ph_lo;
ph_center_rad = (global_ph_hi + global_ph_lo) / 2;
ph_buffer_rad = deg2rad(ph_buffer_deg);
ph_lo_buf     = global_ph_lo + ph_buffer_rad;
ph_hi_buf     = global_ph_hi - ph_buffer_rad;

if global_ph_lo >= global_ph_hi
    error('Phase ranges do not intersect. Check calibration data.');
end
if ph_lo_buf >= ph_hi_buf
    error('ph_buffer_deg (%.1f) exceeds half the common range (%.1f deg).', ...
        ph_buffer_deg, rad2deg(ph_range_rad)/2);
end

fprintf('Phase range (common):    [%+.1f, %+.1f] deg  (%.1f deg span)\n', ...
    rad2deg(global_ph_lo), rad2deg(global_ph_hi), rad2deg(ph_range_rad));
fprintf('Phase range (buffered):  [%+.1f, %+.1f] deg  (%.1f deg span)  centre %+.1f deg\n\n', ...
    rad2deg(ph_lo_buf), rad2deg(ph_hi_buf), ...
    rad2deg(ph_hi_buf - ph_lo_buf), rad2deg(ph_center_rad));

%% =========================
% TARGET DEFINITION
%
% PER-ELEMENT-AWARE AMPLITUDE SCALING
% -----------------------------------
% The goal is to find the single largest global scale k such that no
% element of the normalised matrix exceeds its own calibration ceiling:
%
%   k * |C_norm(r,c)| <= max_amp(r,c)  for all active (r,c)
%
% Rearranging:  k <= max_amp(r,c) / |C_norm(r,c)|  for each active element
%
% Therefore:    k = min over active elements of [ max_amp(r,c) / |C_norm(r,c)| ]
%
% This is tighter than using global_max_amp when some elements have
% larger per-element ceilings, and it is looser (allowing higher overall
% amplitude) when the most restrictive element has a small amplitude
% ratio in the target matrix.
%
% amp_buffer is applied as a final fractional scale-back.  At the
% default value of 0, targets reach the calibration boundary exactly
% and the solver finds the nearest (boundary) calibration point cleanly.
%% =========================

% Element-presence threshold: below this amplitude ratio, an element
% is treated as zero (maximum attenuation) and excluded from the scale
% calculation.  Prevents near-zero elements from driving the scale
% toward infinity.
zero_thresh = 1e-6;

if loadC
    fprintf('Load Compression Matrix\n');
    [file, folder] = uigetfile( ...
        {'*.mat;*.csv;*.txt','Target Files (*.mat, *.csv, *.txt)'}, ...
        'Select Target Compression Matrix');
    target_file = load(fullfile(folder, file));

    C_amp_target = abs(target_file.A);
    C_ph_target  = rad2deg(angle(target_file.A));

    C_norm = target_file.A / max(abs(target_file.A(:)));

    % Phase-aware ceiling: for each element, the achievable amplitude at its
    % specific target phase ? not the global voltage-optimum maximum.
    % Prevents the solver being given a target reachable only at a voltage
    % incompatible with the required phase (amplitude-phase coupling).
    amp_ceil_load = compute_phase_aware_ceiling(Cal_complex, C_norm, N, zero_thresh);
    C_target = scale_to_calibration(C_norm, amp_ceil_load, amp_buffer, zero_thresh);

    check_phase_bounds(C_target, ph_lo_buf, ph_hi_buf, global_ph_lo, global_ph_hi, zero_thresh);

elseif randC
    fprintf('Randomly Generated Compression Matrix\n');

    % Worst-case ceiling: minimum amplitude at i=1 across all calibration
    % voltage points.  Guarantees the amplitude target is achievable at
    % ANY phase, since the solver selects the voltage (and hence phase)
    % independently.  More conservative than phase-aware ceiling but correct
    % for random targets whose phases are unknown at scaling time.
    amp_ceil_rand = compute_worstcase_ceiling(calibration, N);
    wc_global = min(amp_ceil_rand(:));

    rand_amp   = wc_global * rand(N,N);          % amplitude in [0, worst-case ceiling]
    rand_phase = ph_lo_buf + (ph_hi_buf - ph_lo_buf) * rand(N,N);

    C_target     = rand_amp .* exp(1j * rand_phase);
    C_amp_target = abs(C_target);
    C_ph_target  = rad2deg(angle(C_target));

    C_target     = rand_amp .* exp(1j * rand_phase);
    C_amp_target = abs(C_target);
    C_ph_target  = rad2deg(angle(C_target));

elseif useCalVal
    fprintf('Calibration-point Targets (sanity check)\n');

    calix = 40 * ones(N,N);
    caliy =  6 * ones(N,N);

    C_target     = zeros(N,N,'like',1j);
    C_amp_target = zeros(N,N);
    C_ph_target  = zeros(N,N);

    for row = 1:N
        for col = 1:N
            C_target(row,col)     = Cal_complex{row,col}(calix(row,col), caliy(row,col));
            C_amp_target(row,col) = abs(C_target(row,col));
            C_ph_target(row,col)  = rad2deg(angle(C_target(row,col)));
        end
    end

else
    fprintf('Manually Defined Compression Matrix\n');

    % ---- USER-DEFINED COMPRESSION MATRIX --------------------------------
    C_amp_target = eye(N);     % amplitude ratios in [0, 1]
    C_ph_target  = zeros(N);   % phases in degrees; must be within buffered range
    % ---------------------------------------------------------------------

    C_norm   = (C_amp_target ./ max(C_amp_target(:))) ...
             .* exp(1j * deg2rad(C_ph_target));
    amp_ceil_manual = compute_phase_aware_ceiling(Cal_complex, C_norm, N, zero_thresh);
    C_target = scale_to_calibration(C_norm, amp_ceil_manual, amp_buffer, zero_thresh);

    check_phase_bounds(C_target, ph_lo_buf, ph_hi_buf, global_ph_lo, global_ph_hi, zero_thresh);

end

% Display the effective scale relative to global_max_amp so the user
% can see how much of the available dynamic range is being used.
active_scale = max(abs(C_target(:)));
fprintf('Effective target ceiling: %.2f dB  (global max: %.2f dB,  margin: %.2f dB)\n\n', ...
    20*log10(active_scale), 20*log10(global_max_amp), ...
    20*log10(global_max_amp) - 20*log10(active_scale));

%% =========================
% SOLVER OUTPUTS
%% =========================
V                  = zeros(N,N);
Att_dB             = zeros(N,N);
Att_bin            = cell(N,N);
i_indx             = zeros(N,N);
j_indx             = zeros(N,N);
C_complex_hardware = zeros(N,N,'like',1j);

%% =========================
% SOLVE PER ELEMENT
%% =========================
for row = 1:N
    for col = 1:N
        [ii, jj] = find_hw_settings_complex( ...
            Cal_complex{row,col}, C_target(row,col));

        i_indx(row,col) = ii;
        j_indx(row,col) = jj;

        Att_dB(row,col)  = (ii - 1) * att_step;
        Att_bin{row,col} = dec2bin(ii - 1, 6);
        V(row,col)       = (jj - 1) * (V_max / (n_interp - 1));

        C_complex_hardware(row,col) = Cal_complex{row,col}(ii, jj);
    end
end

%% =========================
% ERROR METRICS
%% =========================
amp_error_dB = 20*log10(abs(C_complex_hardware) ./ abs(C_target));
ph_error_deg = rad2deg(angle(C_complex_hardware ./ C_target));
ph_mask      = 20*log10(abs(C_target)) >= ph_report_thresh_dB;
eps_norm     = (C_complex_hardware - C_target) ./ abs(C_target);
RMSE_norm    = sqrt(mean(abs(eps_norm(:)).^2));

%% =========================
% DISPLAY
%% =========================
fprintf('Target compression matrix (mag, dB):\n');    disp(20*log10(abs(C_target)));
fprintf('Hardware compression matrix (mag, dB):\n');  disp(20*log10(abs(C_complex_hardware)));
fprintf('Amplitude error (dB):\n');                   disp(amp_error_dB);
fprintf('  Bias: %+.4f dB    RMSE: %.4f dB\n', ...
    mean(amp_error_dB(:)), sqrt(mean(amp_error_dB(:).^2)));

fprintf('\nTarget compression matrix (phase, deg):\n');   disp(rad2deg(angle(C_target)));
fprintf('Hardware compression matrix (phase, deg):\n');  disp(rad2deg(angle(C_complex_hardware)));
fprintf('Phase error (deg, elements above %d dB):\n', ph_report_thresh_dB);
disp(ph_error_deg .* ph_mask);
if any(ph_mask(:))
    fprintf('  Bias: %+.4f deg    RMSE: %.4f deg\n', ...
        mean(ph_error_deg(ph_mask)), sqrt(mean(ph_error_deg(ph_mask).^2)));
end

fprintf('\nNormalised complex RMSE: %.4f\n', RMSE_norm);
fprintf('\nAttenuator settings (dB):\n');  disp(Att_dB);
fprintf('Phase shifter voltages (V):\n'); disp(V);
fprintf('\nAttenuator binary codes:\n');
for row = 1:N
    fprintf('  ');
    for col = 1:N, fprintf('%-10s', Att_bin{row,col}); end
    fprintf('\n');
end

%% =========================
% ARDUINO EXPORT
%% =========================
if export_to_IDE
    exportToArduinoHeader( ...
        strcat(dothfilename, '.h'), ...
        20*log10(abs(C_target)),       ...
        Att_dB,                        ...
        rad2deg(angle(C_target)),      ...
        V);
end

%% =========================
% LOCAL FUNCTIONS
%% =========================

function C_target = scale_to_calibration(C_norm, amp_ceil, amp_buffer, zero_thresh)
% Scale a normalised complex matrix so the most-constrained active element
% reaches its effective amplitude ceiling (minus optional buffer).
%
% C_norm    : complex [N×N], normalised so max(|C_norm|) = 1
% amp_ceil  : [N×N] effective per-element amplitude ceiling (linear).
%             When called from loadC/manual paths this is the phase-aware
%             ceiling from compute_phase_aware_ceiling().  For randC it is
%             the worst-case (minimum over phase range) ceiling.
% amp_buffer: fractional scale-back in [0,1);  0 = no buffer
% zero_thresh: elements below this ratio treated as zero (no constraint)
%
% Algorithm:
%   k = min over active (r,c) of [ amp_ceil(r,c) / |C_norm(r,c)| ]
%   C_target = C_norm * k * (1 - amp_buffer)

    active = abs(C_norm) > zero_thresh;

    if ~any(active(:))
        error('All elements of C_norm are below zero_thresh. Check input matrix.');
    end

    ratios = amp_ceil(active) ./ abs(C_norm(active));
    k      = min(ratios);

    C_target = C_norm * k * (1 - amp_buffer);

    [r_bind, c_bind] = find(active & (amp_ceil ./ abs(C_norm) <= k * (1 + 1e-9)));
    fprintf('Amplitude scale set by element (%d,%d):  ', r_bind(1), c_bind(1));
    fprintf('target %.2f dB,  ceiling %.2f dB\n', ...
        20*log10(abs(C_target(r_bind(1), c_bind(1)))), ...
        20*log10(amp_ceil(r_bind(1), c_bind(1))));
end

function amp_ceil = compute_phase_aware_ceiling(Cal_complex, C_norm, N, zero_thresh)
% Compute the achievable amplitude ceiling for each element at its specific
% target phase, using the interpolated calibration surface.
%
% The solver operates on the interpolated Cal_complex surface, so the
% ceiling must be computed from the same surface ? not from the original
% coarser calibration data.  Using the original data creates a resolution
% mismatch: the ceiling is computed at the nearest original calibration
% voltage while the solver finds the nearest interpolated voltage, which
% may be at a different operating point with lower amplitude.
%
% For each active element (r,c):
%   1. Determine the target phase from C_norm.
%   2. Scan row 1 (minimum attenuation) of Cal_complex ? the interpolated
%      surface ? to find the voltage index whose phase is nearest to the
%      target phase.
%   3. Use the amplitude at that interpolated point as the ceiling.
%
% This ensures the ceiling reflects what the solver will actually achieve,
% making the scale self-consistent with no residual amplitude error from
% the ceiling computation step.

    global_max_lin = max(cellfun(@(c) max(abs(c(:))), Cal_complex));
    amp_ceil = global_max_lin * ones(N, N);

    for row = 1:N
        for col = 1:N
            if abs(C_norm(row,col)) <= zero_thresh
                continue
            end

            target_ph = angle(C_norm(row,col));

            % Use interpolated Cal_complex row 1 (minimum attenuation).
            % angle() on complex calibration values is consistent with
            % the phase used in the solver's complex-plane distance.
            phases_interp = angle(Cal_complex{row,col}(1,:));   % 1 × n_interp
            [~, nearest_j] = min(abs(phases_interp - target_ph));
            amp_ceil(row,col) = abs(Cal_complex{row,col}(1, nearest_j));
        end
    end
end

function amp_ceil = compute_worstcase_ceiling(calibration, N)
% Amplitude ceiling for randC: minimum amplitude at i=1 across all
% calibration voltage points.  This guarantees the ceiling is achievable
% at ANY phase within the calibration range, regardless of which voltage
% the solver ultimately selects.

    amp_ceil = zeros(N, N);
    for row = 1:N
        for col = 1:N
            raw_amp_row1 = 10.^(calibration(row,col).C_mag(1,:)/20);
            amp_ceil(row,col) = min(raw_amp_row1);
        end
    end
end

function check_phase_bounds(C_target, ph_lo_buf, ph_hi_buf, ...
                             global_ph_lo, global_ph_hi, zero_thresh)
    active        = abs(C_target) > zero_thresh;
    target_phases = angle(C_target);
    out_hard      = active & (target_phases < global_ph_lo | target_phases > global_ph_hi);
    out_buf       = active & ~out_hard & ...
                    (target_phases < ph_lo_buf | target_phases > ph_hi_buf);

    if any(out_hard(:))
        warning(['%d element(s) have target phases OUTSIDE the hardware achievable ' ...
                 'range [%.1f, %.1f] deg and will be clamped by the solver.'], ...
            sum(out_hard(:)), rad2deg(global_ph_lo), rad2deg(global_ph_hi));
    end
    if any(out_buf(:))
        warning(['%d element(s) are within hardware range but outside the symmetric ' ...
                 'buffer zone [%.1f, %.1f] deg.'], ...
            sum(out_buf(:)), rad2deg(ph_lo_buf), rad2deg(ph_hi_buf));
    end
end

function exportToArduinoHeader(filename, A, B, C, D)
    fid = fopen(filename, 'w');
    fprintf(fid, '#ifndef MATRICES_H\n#define MATRICES_H\n\n');
    writeMatrix(fid, 'A', A); writeMatrix(fid, 'B', B);
    writeMatrix(fid, 'C', C); writeMatrix(fid, 'D', D);
    fprintf(fid, '#endif\n');
    fclose(fid);
end

function writeMatrix(fid, name, M)
    [rows, cols] = size(M);
    fprintf(fid, 'float %s[%d][%d] = {\n', name, rows, cols);
    for i = 1:rows
        fprintf(fid, '  {');
        for j = 1:cols
            if j < cols, fprintf(fid, '%.6f, ', M(i,j));
            else,         fprintf(fid, '%.6f',   M(i,j)); end
        end
        if i < rows, fprintf(fid, '},\n'); else, fprintf(fid, '}\n'); end
    end
    fprintf(fid, '};\n\n');
end