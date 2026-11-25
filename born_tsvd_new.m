function rec = born_tsvd_new(E_sc_meas, E_int_doi, dA, freq, current, opts, E_sc_bg)
% BORN_TSVD  Linearized Born / DBIM step using Truncated SVD.
%   If E_sc_bg is provided (scattered field from current DBIM background),
%   inverts the data residual (E_sc_meas - E_sc_bg) using distorted fields
%   E_int_doi to compute the contrast update delta_eps.

if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'threshold'), opts.threshold = 0.9999; end
if ~isfield(opts, 'plotSVD'),   opts.plotSVD   = false;  end
if ~isfield(opts, 'alpha'),     opts.alpha     = 1.0;    end
if nargin < 7 || isempty(E_sc_bg)
    E_sc_bg = zeros(size(E_sc_meas), 'like', E_sc_meas);
end

eps0  = 8.854187817e-12;
omega = 2*pi*freq;

%% Data residual vector (reciprocal-unique pairs n <= m)
E_res  = E_sc_meas - E_sc_bg;
N_ant  = size(E_res, 1);
[n, m] = find(triu(true(N_ant)));          % n = transmitter, m = receiver
e      = E_res(sub2ind([N_ant, N_ant], m, n));

%% Distorted data operator K
E_tx = E_int_doi(:, n);                    % Ncell-by-Npairs
E_rx = E_int_doi(:, m);
K    = (1j*omega*eps0/current) * (E_tx.' .* E_rx.') .* dA(:).';

%% Truncated SVD (recomputed each DBIM iteration)
[U, S, V] = svd(K, 'econ');
s          = diag(S);
cum_energy = cumsum(s.^2) / max(sum(s.^2), eps);
Nt         = find(cum_energy >= opts.threshold, 1, 'first');
if isempty(Nt), Nt = numel(s); end

if opts.plotSVD
    plot_svd_analysis(S);
end

x = opts.alpha * (V(:,1:Nt) * ((U(:,1:Nt)' * e) ./ s(1:Nt)));

rec.delta_eps = x;
rec.Nt        = Nt;
rec.s         = s;
rec.K         = K;
rec.res_norm  = norm(e);
end