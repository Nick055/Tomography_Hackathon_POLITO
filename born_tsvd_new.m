function rec = born_tsvd_new(E_sc_meas, E_int_doi, dA, freq, current, opts, E_sc_bg)
% BORN_TSVD  Linearized Born / DBIM step using Truncated SVD.
%   If E_sc_bg is provided (scattered field from current DBIM background),
%   inverts the data residual (E_sc_meas - E_sc_bg) using distorted fields
%   E_int_doi to compute the contrast update delta_eps.

if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'regMethod'), opts.regMethod = 'tsvd'; end   % 'tsvd' or 'tikhonov'
if ~isfield(opts, 'threshold'), opts.threshold = 0.999; end
if ~isfield(opts, 'lambda'),    opts.lambda    = 1e-2;   end   % Tikhonov parameter
if ~isfield(opts, 'relLambda'), opts.relLambda = true;   end   % Scale lambda by s(1)
if ~isfield(opts, 'plotSVD'),   opts.plotSVD   = false;  end
if ~isfield(opts, 'alpha'),     opts.alpha     = 1.0;    end
if nargin < 7 || isempty(E_sc_bg)
    E_sc_bg = zeros(size(E_sc_meas), 'like', E_sc_meas);
end

eps0  = 8.854187817e-12;
omega = 2*pi*freq;

%% Data residual vector (measured pairs, each reciprocal pair once: m >= n)
E_res  = E_sc_meas - E_sc_bg;
N_ant  = size(E_res, 1);

if isfield(opts, 'mask') && ~isempty(opts.mask)
    sel = tril(logical(opts.mask));      % measured pairs only (uses the mask)
else
    sel = tril(true(N_ant));             % all pairs
end
[m, n] = find(sel);                      % m = receiver (row), n = transmitter (col), m >= n

e      = E_res(sub2ind([N_ant, N_ant], m, n));

%% Distorted data operator K
E_tx = E_int_doi(:, n);                    % Ncell-by-Npairs
E_rx = E_int_doi(:, m);
K    = (1j*omega*eps0/current) * (E_tx.' .* E_rx.') .* dA(:).';

%% Truncated SVD (recomputed each DBIM iteration)
%% Regularized inversion via SVD (recomputed each DBIM iteration)
[U, S, V] = svd(K, 'econ');
s = diag(S);

if opts.plotSVD
    plot_svd_analysis(S);
end

switch lower(opts.regMethod)
    case 'tsvd'
        cum_energy = cumsum(s.^2) / max(sum(s.^2), eps);
        Nt         = find(cum_energy >= opts.threshold, 1, 'first');
        if isempty(Nt), Nt = numel(s); end

        x = opts.alpha * (V(:,1:Nt) * ((U(:,1:Nt)' * e) ./ s(1:Nt)));
        rec.Nt = Nt;

    case 'tikhonov'
        if opts.relLambda
            lam = opts.lambda * s(1);   % Relative to largest singular value
        else
            lam = opts.lambda;          % Absolute lambda
        end

        filt = s ./ (s.^2 + lam^2);
        x    = opts.alpha * (V * (filt .* (U' * e)));
        rec.lambda = lam;

    otherwise
        error('Unknown opts.regMethod: %s. Use ''tsvd'' or ''tikhonov''.', opts.regMethod);
end

rec.delta_eps = x;
rec.s         = s;
rec.K         = K;
rec.res_norm  = norm(e);
end