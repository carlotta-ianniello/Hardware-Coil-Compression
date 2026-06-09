clear

%save_settings_txt = 0; % save settings in text file
export_to_IDE = 1;% Generate .h file to export into Arduino IDE
%dothfilenamepath = 'C:\Users\iannic01\Personal\Research\Arduino_test\modular_board_test\4x4_matrix\4x4_calibration_mod\calibration_mod\4x4_validation\4x4_validation_mod\';
dothfilenamepath = 'C:\Users\iannic01\Personal\Research\Arduino_test\modular_board_test\4x4_matrix\4x4_matrix_settings';
dothfilename = [dothfilenamepath 'setting_compression_matrix_ch1357'];

%% =========================
% USER TARGET (4x4 INPUT)
%% =========================
loadC = 1;    % Load Calibration matrix
randC = 0;    % Generate random Calibration matrix with Att = [-40 -16] and Ph = [0 180] (used for validation)
useCalVal = 0; % Set as target Calibration matrix a combination of settings that was measured during calibration (sanity check)

%if (load_C && rand_C && useCalVal) || (load_C && rand_C) || (rand_C && useCalVal) || (load_C && useCalVal)
if loadC+randC+useCalVal>1
    error('Please select either loadC, randC or useCalVal');
end


if (loadC)
    fprintf('Load Compression Matrix\n');
    [file, folder] = uigetfile( ...
        {'*.mat;*.csv;*.txt','Target Files (*.mat, *.csv, *.txt)'}, ...
        'Select Target Compression Matrix');
    target_file = load(fullfile(folder, file));

    C_amp_target = abs(target_file.A);
    C_ph_target = rad2deg(angle(target_file.A));

    C_amp_rel = C_amp_target./max(C_amp_target(:));
    C_ph_rel  = C_ph_target;

elseif randC
    fprintf('Randomly Generated Compression Matrix\n');
    % range_amp = [0 1];
    % C_amp_target = range_amp(1) + (range_amp(2) - range_amp(1)) * rand(4,4);
    range_amp = [16 40];
    C_amp_target = range_amp(1) + (range_amp(2) - range_amp(1)) * rand(4,4);
    range_ph = [0 180];
    C_ph_target = range_ph(1) + (range_ph(2) - range_ph(1)) * rand(4,4);

    C_amp_rel = - C_amp_target;
    C_ph_rel  = C_ph_target;

elseif useCalVal
    fprintf('Compression Matrix with combinations from Calibration\n');
    C_amp_rel = [];
    C_ph_rel = [];
else

    fprintf('Manually defined Compression Matrix\n');

    % -- USER DEFINED COMPRESSION MATRIX GOES HERE!--------
    C_amp_target = eye(4);
    C_ph_target  = ones(4);
    %------------------------------------------------------

    C_amp_rel = C_amp_target./max(C_amp_target(:));
    C_ph_rel  = C_ph_target;
end



%% =========================
% PATH + LOAD CALIBRATION
%% =========================
path = 'C:\Users\iannic01\Personal\Research\Arduino_test\modular_board_test\4x4_matrix\4x4_calibration_mod\calibration_mod\';
addpath(genpath(path));

load([path 'calibration_4x4.mat']);  % structured calibration

n_target = 64;

Cal_mag = cell(4,4);
Cal_ph  = cell(4,4);

%% =========================
% EXTRACT + NORMALIZE CALIBRATION TO 64x64
%% =========================
for row = 1:4
    for col = 1:4

        mag = calibration(row,col).C_mag - 1.12;    % added bias prom prev meas
        ph  = calibration(row,col).C_ph -5.4;

        [m, n] = size(mag);

        x_old = linspace(1, n, n);
        x_new = linspace(1, n, n_target);

        y_old = linspace(1, m, m);
        y_new = linspace(1, m, n_target);

        % unwrap phase
        ph = mod(ph + 180, 360) - 180;

        % interpolate columns
        mag_i = interp1(x_old, mag.', x_new, 'linear', 'extrap').';
        ph_i  = interp1(x_old, ph.',  x_new, 'linear', 'extrap').';

        % interpolate rows
        Cal_mag{row,col} = interp1(y_old, mag_i, y_new, 'linear', 'extrap');
        Cal_ph{row,col}  = interp1(y_old, ph_i,  y_new, 'linear', 'extrap');

        Cal_ph{row,col} = mod(Cal_ph{row,col} + 180, 360) - 180;
    end
end

%% =========================
% GLOBAL BOUNDARIES
%% =========================
minAtt = inf(4,4);
maxAtt = inf(4,4);
minPh  = inf(4,4);
maxPh  = inf(4,4);

for row = 1:4
    for col = 1:4
        minAtt(row,col) = max(Cal_mag{row,col}(:));
        maxAtt(row,col) = min(Cal_mag{row,col}(:));

        minPh(row,col) = min(Cal_ph{row,col}(:));
        maxPh(row,col) = max(Cal_ph{row,col}(:));
    end
end

boundariesAtt = [min(minAtt(:)) max(maxAtt(:))];
boundariesPh  = [max(minPh(:)) min(maxPh(:))];

%% =========================
% SHIFT PHASE TO POSITIVE RANGE
%% =========================
phase_offset = abs(min(minPh(:)));

for row = 1:4
    for col = 1:4
        Cal_ph{row,col} = Cal_ph{row,col} + phase_offset;
    end
end

%% =========================
% TARGET FORMATTING
%% =========================
if ~(randC)
    %C_amp = -(-ceil(boundariesAtt(1)) + 1 - db(C_amp_rel,'power'));
    %C_amp = -(-ceil(boundariesAtt(1)) + 1 - db(C_amp_rel,'voltage'));
    C_amp = -(-ceil(boundariesAtt(1)) - db(C_amp_rel,'voltage'));

elseif randC
    C_amp = C_amp_rel;
end
%C_ph  = C_ph_rel + 10 + phase_offset;
C_ph  = C_ph_rel + 10 - min(C_ph_rel(:));

C_amp(isinf(C_amp)) = -45;



% if randC
%     C_amp = C_amp_rel;
%     C_ph  = C_ph_rel +10 - min(C_ph_rel(:));
% elseif useCalVal
%     for xx = 1:4
%         for yy = 1:4
%             calix = 40;%randi([1 64],1,1);
%             caliy = 6 ;%randi([1 7],1,1);
%             C_amp(xx,yy) = calibration(xx,yy).C_mag(calix,caliy);
%             C_ph(xx,yy) = calibration(xx,yy).C_ph(calix,caliy);
%             C_ph  = C_ph +10 - min(C_ph(:));
%         end 
%     end 
% else
%     C_amp = -(-ceil(boundariesAtt(1)) + 1 - db(C_amp_rel,'power'));
%     C_ph  = C_ph_rel +10 - min(C_ph_rel(:));
% end
% 
% C_amp(isinf(C_amp)) = -45;

%% =========================
% SOLVER OUTPUTS
%% =========================
V = zeros(4,4);
Att_dB = zeros(4,4);
Att_bin = cell(4,4);

i_indx = zeros(4,4);
j_indx = zeros(4,4);

Att_dB_hardware = zeros(4,4);
Ph_deg_hardware = zeros(4,4);

%% =========================
% SOLVE PER CHANNEL
%% =========================
for row = 1:4
    for col = 1:4

        target_mag_dB = C_amp(row,col);
        target_phase_deg = C_ph(row,col);

        [ii, jj, V_temp, Att_dB_temp, Att_bin_temp] = ...
            find_hw_settings_dB_phase( ...
                Cal_mag{row,col}, ...
                Cal_ph{row,col}, ...
                target_mag_dB, ...
                target_phase_deg);

        i_indx(row,col) = ii;
        j_indx(row,col) = jj;

        V(row,col) = V_temp;
        Att_dB(row,col) = Att_dB_temp;
        Att_bin{row,col} = Att_bin_temp;

        Att_dB_hardware(row,col) = Cal_mag{row,col}(ii,jj);
        Ph_deg_hardware(row,col) = Cal_ph{row,col}(ii,jj);

    end
end

%% =========================
% RESULTS
%% =========================
disp('Target amplitude matrix:');
disp(C_amp);

disp('Hardware amplitude matrix:');
disp(Att_dB_hardware);

disp('Attenuator settings (dB):');
disp(Att_dB);

disp('Attenuator binary settings:');
disp(Att_bin);

disp('Target phase matrix:');
disp(C_ph);

disp('Hardware phase matrix:');
disp(Ph_deg_hardware);

disp('Phase shifter voltages:');
disp(V);

%% =========================
% EXPORT OPTIONS
%% =========================

if export_to_IDE
    %t = datetime("now");
    %timestamp = datestr(t, "yyyy-mm-dd_HH-MM-SS");
    %filename = (path + "4x4_validation\settings_" + timestamp +".h");
    %mkdir C:\Users\iannic01\Personal\Research\Arduino_test\modular_board_test\4x4_matrix\4x4_calibration_mod\calibration_mod\4x4_validation\mix_mat_validation ...
    %setting11
    filename = (dothfilename + ".h");
    %filename = "settings_" + timestamp +".h";
    exportToArduinoHeader(filename, C_amp, Att_dB, C_ph, V);
end


%%
    function exportToArduinoHeader(filename, A, B, C, D)
    fid = fopen(filename, 'w');

    fprintf(fid, '#ifndef MATRICES_H\n#define MATRICES_H\n\n');

    writeMatrix(fid, 'A', A);
    writeMatrix(fid, 'B', B);
    writeMatrix(fid, 'C', C);
    writeMatrix(fid, 'D', D);

    fprintf(fid, '#endif\n');
    fclose(fid);
end

function writeMatrix(fid, name, M)
    [rows, cols] = size(M);
    fprintf(fid, 'float %s[%d][%d] = {\n', name, rows, cols);

    for i = 1:rows
        fprintf(fid, '  {');
        for j = 1:cols
            if j < cols
                fprintf(fid, '%.6f, ', M(i,j));
            else
                fprintf(fid, '%.6f', M(i,j));
            end
        end
        if i < rows
            fprintf(fid, '},\n');
        else
            fprintf(fid, '}\n');
        end
    end

    fprintf(fid, '};\n\n');
end