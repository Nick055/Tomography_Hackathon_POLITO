% Main script
clc;
close all;
clear;

resultsDir = fullfile("Examples/Circular_scattered_field/results");

operatorName  = "operator";   % E_inc on the reconstruction mesh
stateName     = "total";      % measured state (with target)
referenceName = "incident";   % reference state (no target)

operator  = load_fields(operatorName,  resultsDir);
state     = load_fields(stateName,     resultsDir);
reference = load_fields(referenceName, resultsDir);

[E_sc, truth] = differential_data(state, reference, operator);

doi = operator.mesh.cellsDOI;

opts.threshold = 0.999;
opts.plotSVD   = true;
rec = born_tsvd(E_sc, operator.Ecells(doi, :), operator.mesh.area(doi), ...
                operator.freq, operator.current, opts);

verts = operator.mesh.verts;
cells = operator.mesh.cells(doi, :);

plotContrast(verts, cells, real(rec.delta_eps), "Recovered $\mathrm{Re}\{\Delta\varepsilon_r\}$");
plotContrast(verts, cells, imag(rec.delta_eps), "Recovered $\mathrm{Im}\{\Delta\varepsilon_r\}$");

img_norm = abs(rec.delta_eps) ./ max(abs(rec.delta_eps));
plotContrast(verts, cells, img_norm, "Normalised $|\Delta\varepsilon_r|$");

binThreshold = 10^(-5/20);
plotContrast(verts, cells, double(img_norm > binThreshold), "Binarised $|\Delta\varepsilon_r|$ ($-5$ dB)");

if ~isempty(truth)
    plotContrast(truth.verts, truth.cells, real(truth.contrast), "True $\mathrm{Re}\{\Delta\varepsilon_r\}$");
    plotContrast(truth.verts, truth.cells, imag(truth.contrast), "True $\mathrm{Im}\{\Delta\varepsilon_r\}$");
end