function [i_best, j_best, voltage, attenuation_dB, attenuator_bin] = ...
    find_hw_settings_dB_phase(C_mag, C_ph, target_mag_dB, target_phase_deg)

    % % Check matrix sizes
    % if ~isequal(size(C_mag), [64, 13]) || ~isequal(size(C_ph), [64, 13])
    %     error('C_mag and C_ph must be 64x13 matrices.');
    % end

    [n,m] = size(C_mag);

    % --- Convert C_mag (negative dB) to linear magnitude
    mag_lin = 10.^(C_mag / 20);

    % --- Convert C_ph (degrees) to radians
    phase_rad = deg2rad(C_ph);

    % --- Build complex lookup table
    C = mag_lin .* exp(1j * phase_rad);

    % --- Convert target mag_dB to linear magnitude
    target_mag_lin = 10^(target_mag_dB / 20);

    % --- Convert target phase (deg) to radians
    target_phase_rad = deg2rad(target_phase_deg);

    % --- Build target complex number
    target_complex = target_mag_lin * exp(1j * target_phase_rad);

    % % --- Find closest entry in complex plane
    % diffMatrix = C - target_complex;
    % distMatrix = abs(diffMatrix);
    %
    % [~, minIdx] = min(distMatrix(:));
    % [i_best, j_best] = ind2sub(size(C), minIdx);
    % Magnitude error (in dB)
    mag_error = abs(C_mag - target_mag_dB);

    % Phase error (wrapped properly!)
    phase_error = abs(C_ph - target_phase_deg);
    phase_error = min(phase_error, 360 - phase_error); % wrap-around fix

    % Normalize 
    mag_weight = 1;
    phase_weight = 1/5; % tune this!

    % Combined cost
    costMatrix = mag_weight * mag_error + phase_weight * phase_error;

    % Find minimum
    [~, minIdx] = min(costMatrix(:));
    [i_best, j_best] = ind2sub(size(C_mag), minIdx);

    % --- Convert j index to voltage (0–12 V, m steps)
    voltage = (j_best - 1) * (12 / (m - 1));

    % --- Convert i index to attenuation in dB
    attenuation_dB = (i_best - 1) * 0.5;

    % --- Binary code for attenuator
    attenuator_bin = dec2bin(i_best - 1, 6);

    % --- Display results
    % fprintf('Best match: i = %d, j = %d\n', i_best, j_best);
    % fprintf('Phase Shifter Voltage: %.2f V\n', voltage);
    % fprintf('Attenuation: %.1f dB\n', attenuation_dB);
    % fprintf('Attenuator binary code: %s\n', attenuator_bin);
end
