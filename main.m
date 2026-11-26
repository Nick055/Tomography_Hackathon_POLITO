% Main script
clc;
close all;
clear;

setup_osmi;
resultsDir = fullfile("Datasets/default_12antennas");
maxIter = 15;
min_err = 1e-2;
include_reflections = true;
solver = "DBIM"; %BIM or DBIM

%% LOAD DATA

operatorName  = "operator";   % E_inc on the reconstruction mesh
stateName     = "total";      % measured state (with target)
referenceName = "incident";   % reference state (no target)

operator  = load_fields(operatorName,  resultsDir);
state     = load_fields(stateName,     resultsDir);
reference = load_fields(referenceName, resultsDir);

%% SETUP MASK: which Tx-Rx pairs are measured
N    = size(operator.Eprobes, 1);        % number of antennas
idx  = (1:N).';
mask = mod(idx,2) ~= mod(idx.',2);       % odd-even pairs (not all antennas pairs were measured)
if include_reflections
    mask(1:N+1:end) = true;              % diagonal = reflections (monostatic)
end
opts.mask = mask;                        % this will be passed into "born_tsvd_new.m"

%% COMPUTE BIM
[E_sc, truth] = differential_data(state, reference, operator);

doi   = operator.mesh.cellsDOI;
areas = operator.mesh.area(doi);
E_inc_doi = operator.Ecells(doi,:);

opts.threshold = 0.99;
opts.plotSVD   = false;

binThreshold = 10^(-5/20);

if solver == "BIM"
    
    E_int = E_inc_doi;
    delta_eps_prev = zeros(numel(doi),1);
    
    for ii = 1:maxIter
        % 1. Linearized inverse step (TSVD) using current internal field E_int
        rec = born_tsvd_new(E_sc, E_int, areas, operator.freq, operator.current, opts);
        delta_eps = rec.delta_eps;
    
        % 2. Check parameter convergence
        % TODO: think of better metric for relChange
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
   
    plotContrast(verts, cells, double(img_norm > binThreshold), "Binarised $|\Delta\varepsilon_r|$ ($-5$ dB)");
    
    if ~isempty(truth)
        plotContrast(truth.verts, truth.cells, real(truth.contrast), "True $\mathrm{Re}\{\Delta\varepsilon_r\}$");
        plotContrast(truth.verts, truth.cells, imag(truth.contrast), "True $\mathrm{Im}\{\Delta\varepsilon_r\}$");
    end

%% COMPUTE DBIM
elseif solver == "DBIM"
    alpha          = 0.7;          % Relaxation / damping factor in (0, 1]
    perm_bg_doi    = operator.perm(doi);
    perm_bg_doi    = perm_bg_doi(:);
    perm_dbim      = operator.perm(:);
    E_int_dbim     = E_inc_doi;
    E_sc_bg        = zeros(size(E_sc), 'like', E_sc);

    best_res_norm  = inf;
    best_perm_dbim = perm_dbim;

    opts.regMethod = 'tikhonov';
    opts.relLambda = true;
    opts.lambda = 1e-1;

    for ii = 1:maxIter
        % 1. Linearized DBIM step using current distorted internal field and data residual
        rec_dbim = born_tsvd_new(E_sc, E_int_dbim, areas, operator.freq, operator.current, opts, E_sc_bg);
        delta_eps_step = rec_dbim.delta_eps;

        % Track lowest-residual iterate
        if isfield(rec_dbim, 'res_norm') && rec_dbim.res_norm < best_res_norm
            best_res_norm  = rec_dbim.res_norm;
            best_perm_dbim = perm_dbim;
        end

        % 2. Update cumulative permittivity in DOI (damped) and check profile convergence
        perm_prev_doi  = perm_dbim(doi);
        perm_dbim(doi) = perm_dbim(doi) + alpha * delta_eps_step(:);

        relChange_dbim = norm(perm_dbim(doi) - perm_prev_doi) / max(norm(perm_dbim(doi)), eps);
        if isfield(rec_dbim, 'res_norm')
            fprintf("[DBIM] Iter %d: relChange = %.9e | res_norm = %.9e\n", ...
                ii, relChange_dbim, rec_dbim.res_norm);
        else
            fprintf("[DBIM] Iter %d: relChange = %.9e\n", ii, relChange_dbim);
        end

        if ii > 1 && relChange_dbim < min_err
            break;
        end

        % 3. Update distorted internal field and forward-modeled antenna scattered field
        sol_dbim = solve_fields(operator.mesh.verts, ...
            operator.mesh.cells, ...
            perm_dbim, ...
            operator.freq, ...
            operator.mesh.probeNodes, ...
            operator.eps_bg, ...
            operator.current, ...
            operator.mesh.cellsDOI);

        E_int_dbim = sol_dbim.Ecells(doi, :);

        if isfield(sol_dbim, 'E_sc_probes')
            E_sc_bg = sol_dbim.E_sc_probes .* mask;
        elseif isfield(sol_dbim, 'Eprobes')
            E_sc_bg = (sol_dbim.Eprobes - reference.Eprobes) .* mask;
        else
            error('Specify how to extract antenna scattered field E_sc_bg from sol_dbim.');
        end
    end

    % Evaluate final iterate against lowest-residual iterate
    if isfinite(best_res_norm)
        sol_check = solve_fields(operator.mesh.verts, ...
            operator.mesh.cells, ...
            perm_dbim, ...
            operator.freq, ...
            operator.mesh.probeNodes, ...
            operator.eps_bg, ...
            operator.current, ...
            operator.mesh.cellsDOI);
        if isfield(sol_check, 'E_sc_probes')
            final_res = norm(tril(mask) .* (E_sc - sol_check.E_sc_probes), 'fro');
        elseif isfield(sol_check, 'Eprobes')
            final_res = norm(tril(mask) .* (E_sc - (sol_check.Eprobes - reference.Eprobes)), 'fro');
        else
            final_res = inf;
        end
        if best_res_norm < final_res
            fprintf("[DBIM] Restoring best iterate (res_norm %.9e < %.9e)\n", best_res_norm, final_res);
            perm_dbim = best_perm_dbim;
        end
    end

    verts = operator.mesh.verts;
    cells = operator.mesh.cells(doi,:);
    delta_eps_dbim = perm_dbim(doi) - perm_bg_doi;

    plotContrast(verts, cells, real(delta_eps_dbim), "DBIM: Recovered $\mathrm{Re}\{\Delta\varepsilon_r\}$");
    plotContrast(verts, cells, imag(delta_eps_dbim), "DBIM: Recovered $\mathrm{Im}\{\Delta\varepsilon_r\}$");

    img_norm_dbim = abs(delta_eps_dbim) ./ max(abs(delta_eps_dbim));
    plotContrast(verts, cells, img_norm_dbim, "DBIM: Normalised $|\Delta\varepsilon_r|$");
    plotContrast(verts, cells, double(img_norm_dbim > binThreshold), "DBIM: Binarised $|\Delta\varepsilon_r|$ ($-5$ dB)");
else
    error("Invalid solver.")
end


%% Matrix structure 
% entries

% N = 12;
% M = mod(1:N,2) ~= mod((1:N).',2); % Mask

% imagesc(M); axis image; colormap(gray); colorbar; xlabel('Rx'); ylabel('Tx'); % plot the mask to check it
% xticks(1:size(M,1)); yticks(1:size(M,1)); % plot all the thicks on the rulers
% K= magic(12);
%Test = K .*M;% logical multiplication 