clc;
close all;
clear all;

setup_osmi;

simFile = fullfile("Datasets/Measurements_lab/domain.mat");
expFile = fullfile("Datasets/Measurements_lab/meas_12_pla_circle_pos_1.s12p");
expBgdFile = fullfile("Datasets/Measurements_lab/meas_11_background.s12p");
maxIter = 15;
min_err = 1e-2;
binThreshold = 10^(-5/20);
include_reflections = false; %include diagonal of S matrix (S11, S22 etc)
include_both_triangles = false;
solver = "DBIM"; %BIM or DBIM

opts.regMethod = "tikhonov"; %'tsvd' or 'tikhonov'
opts.relLambda = true; %scale lambda by largest singular value s(1)
opts.lambda = 0.1;
opts.includeReflections = include_reflections; 
opts.plotSVD = false;
opts.threshold = 0.99;
alpha = 0.5;


%% LOAD DATA

simMat = load(simFile);
sim_E = simMat.fields.Eprobes;

expSP = sparameters(expFile);
spFreqs = expSP.Frequencies;
[~,idx] = min(abs(spFreqs - simMat.fields.freq));
exp_S = expSP.Parameters(:,:,idx);

expBgdSP = sparameters(expBgdFile);
expBgd_S = expBgdSP.Parameters(:,:,idx);

%% CALIBRATE EXPERIMENTALLY-OBTAINED S-PARAMETERS

calib = load("calibration.mat");
C = calib.C;
C(isnan(C)) = 0;
exp_E = C .* (exp_S - expBgd_S);


%% SETUP MASK FOR S PARAMETERS
N = size(exp_S, 1);
[r, c] = ndgrid(1:N, 1:N);
mask = ones(N, N);
mask(mod(r, 2) == mod(c, 2)) = 0;
if include_reflections
    mask(logical(eye(N))) = 1;
end
if ~include_both_triangles
    mask(r < c) = 0;
end

figure;
imagesc(mask); axis image; colormap(gray); colorbar;
xlabel('Rx'); ylabel('Tx');
xticks(1:N); yticks(1:N);
title('Active Tx/Rx Pair Mask');

% "Clutter Removal"
E_sc = exp_E .* mask;
truth = [];


%% DBIM SOLVER
doi = simMat.fields.mesh.cellsDOI;
areas = simMat.fields.mesh.area(doi);
E_inc_doi = simMat.fields.Ecells(doi,:);
if solver == "DBIM"    
    perm_bg_doi = simMat.fields.perm(doi);
    perm_bg_doi = perm_bg_doi(:);
    perm_dbim = simMat.fields.perm(:);
    E_int_dbim = E_inc_doi;
    E_sc_bg = zeros(size(E_sc), 'like', E_sc);

    best_res_norm = inf;
    best_perm_dbim = perm_dbim;

    for ii = 1:maxIter
        % 1. Linearized DBIM step using current distorted internal field and data residual
        rec_dbim = born_tsvd_new(E_sc, ...
            E_int_dbim, areas, simMat.fields.freq, ...
            simMat.fields.current, opts, E_sc_bg);
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
        sol_dbim = solve_fields(simMat.fields.mesh.verts, ...
            simMat.fields.mesh.cells, ...
            perm_dbim, ...
            simMat.fields.freq, ...
            simMat.fields.mesh.probeNodes, ...
            simMat.fields.eps_bg, ...
            simMat.fields.current, ...
            simMat.fields.mesh.cellsDOI);

        E_int_dbim = sol_dbim.Ecells(doi, :);

        if isfield(sol_dbim, 'E_sc_probes')
            E_sc_bg = sol_dbim.E_sc_probes .* mask;
        elseif isfield(sol_dbim, 'Eprobes')
            E_sc_bg = (sol_dbim.Eprobes - simMat.fields.Eprobes) .* mask;
        else
            error('Specify how to extract antenna scattered field E_sc_bg from sol_dbim.');
        end
    end

    % Evaluate final iterate against lowest-residual iterate
    if isfinite(best_res_norm)
        sol_check = solve_fields(simMat.fields.mesh.verts, ...
            simMat.fields.mesh.cells, ...
            perm_dbim, ...
            simMat.fields.freq, ...
            simMat.fields.mesh.probeNodes, ...
            simMat.fields.eps_bg, ...
            simMat.fields.current, ...
            simMat.fields.mesh.cellsDOI);
        if isfield(sol_check, 'E_sc_probes')
            final_res = norm((E_sc - sol_check.E_sc_probes .* mask), 'fro');
        elseif isfield(sol_check, 'Eprobes')
            final_res = norm((E_sc - (sol_check.Eprobes - simMat.fields.Eprobes) .* mask), 'fro');
        else
            final_res = inf;
        end
        if best_res_norm < final_res
            fprintf("[DBIM] Restoring best iterate (res_norm %.9e < %.9e)\n", best_res_norm, final_res);
            perm_dbim = best_perm_dbim;
        end
    end

    verts = simMat.fields.mesh.verts;
    cells = simMat.fields.mesh.cells(doi, :);
    delta_eps_dbim = perm_dbim(doi) - perm_bg_doi;

    plotContrast(verts, cells, real(delta_eps_dbim), "DBIM: Recovered $\mathrm{Re}\{\Delta\varepsilon_r\}$");
    plotContrast(verts, cells, imag(delta_eps_dbim), "DBIM: Recovered $\mathrm{Im}\{\Delta\varepsilon_r\}$");

    img_norm_dbim = abs(delta_eps_dbim) ./ max(abs(delta_eps_dbim));
    plotContrast(verts, cells, img_norm_dbim, "DBIM: Normalised $|\Delta\varepsilon_r|$");
    plotContrast(verts, cells, double(img_norm_dbim > binThreshold), "DBIM: Binarised $|\Delta\varepsilon_r|$ ($-5$ dB)");

    if ~isempty(truth)
        plotContrast(truth.verts, truth.cells, real(truth.contrast), "True $\mathrm{Re}\{\Delta\varepsilon_r\}$");
        plotContrast(truth.verts, truth.cells, imag(truth.contrast), "True $\mathrm{Im}\{\Delta\varepsilon_r\}$");
    end
end

