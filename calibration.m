clc;
close all;
clear;
setup_osmi;

simBgdFile = fullfile("Datasets/Measurements_lab/domain.mat");
expBgdFile = fullfile("Datasets/Measurements_lab/meas_1_background.s12p");
outputFile = fullfile("calibration.mat");
expFreq = 1e9;

simBgd = load(simBgdFile);
simBgd_E = simBgd.fields.Eprobes;

expBgdSP = sparameters(expBgdFile);
spFreqs = expBgdSP.Frequencies;
[~,idx] = min(abs(spFreqs - expFreq));
expBgd_S = expBgdSP.Parameters(:,:,idx);

% Compute calibration factor and mark unmeasured / invalid ports as NaN
C = simBgd_E ./ expBgd_S;
invalidMask = ~isfinite(C) | isnan(C) | (expBgd_S == 0);
C(invalidMask) = NaN;

validMask = ~isnan(C);
save(outputFile,"C","expFreq","validMask");

magC = abs(C);

figure;
ax = axes;
hImg = imagesc(ax, magC);

% Make NaN entries transparent so the white axes background shows through
set(hImg, 'AlphaData', ~isnan(magC));
set(ax, 'Color', 'w');
colorbar;
xlabel("Receiver Port");
ylabel("Transmitter Port");
title("abs(C) with Unmeasured Switch-Matrix Ports in White");