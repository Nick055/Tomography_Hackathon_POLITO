% Main script
clc;
close all;
clear;

resultsDir = fullfile("Examples/Circular_scattered_field/results");
maxIter = 15;
min_err = 1e-2;

operatorName  = "operator";   % E_inc on the reconstruction mesh
stateName     = "total";      % measured state (with target)
referenceName = "incident";   % reference state (no target)

operator  = load_fields(operatorName,  resultsDir);
state     = load_fields(stateName,     resultsDir);
reference = load_fields(referenceName, resultsDir);

[E_sc, truth] = differential_data(state, reference, operator);

doi   = operator.mesh.cellsDOI;
areas = operator.mesh.area(doi);
E_inc_doi = operator.Ecells(doi,:);

opts.threshold = 0.999;
opts.plotSVD   = false;

E_int = E_inc_doi;
delta_eps_prev = zeros(numel(doi),1);

for ii = 1:maxIter
    % 1. Linearized inverse step (TSVD) using current internal field E_int
    rec = born_tsvd(E_sc, E_int, areas, operator.freq, operator.current, opts);
    delta_eps = rec.delta_eps;

    % 2. Check parameter convergence
    relChange = norm(delta_eps - delta_eps_prev) / max(norm(delta_eps), eps);
    fprintf("Iter %d: relChange = %.9e\n", ii, relChange);
    if ii > 1 && relChange < min_err
        break;
    end
    delta_eps_prev = delta_eps;
    perm_full = operator.perm(:);
    perm_full(doi) = perm_full(doi) + delta_eps(:);

    % 3. Update internal total field inside DOI for next iteration
    sol = solve_fields(operator.mesh.verts, ...
        operator.mesh.cells, ...
        perm_full, ...
        operator.freq, ...
        operator.mesh.probeNodes, ...
        operator.eps_bg, ...
        operator.current, ...
        operator.mesh.cellsDOI);
    E_int = sol.Ecells(doi,:);
end

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


%% Matrix structure 
% entries

N = 12;
M = mod(1:N,2) ~= mod((1:N).',2); % Mask

imagesc(M); axis image; colormap(gray); colorbar; xlabel('Rx'); ylabel('Tx'); % plot the mask to check it
xticks(1:size(M,1)); yticks(1:size(M,1)); % plot all the thicks on the rulers
% K= magic(12);
Test = K .*M;% logical multiplication 