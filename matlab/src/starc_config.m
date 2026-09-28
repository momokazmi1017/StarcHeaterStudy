function cfg = starc_config(p, coat, loadCase, A)
%STARC_CONFIG  Build the N-by-1 case arrays used by STARC_SIMULATE.
%
%   cfg = starc_config(p, coat, loadCase, A)
%
%   coat     : struct with fields alpha, eps_c, eps_h, Tt_C, w
%              (a single coating from p.coat, or column vectors, e.g. from
%              STARC_STACK_COATINGS or a Monte Carlo draw)
%   loadCase : struct with Qint_sun, Qint_ecl [W] (e.g. p.hot or p.cold)
%   A        : radiator area [m^2]
%
%   Every input may be scalar or length N; scalars are expanded. The node's
%   thermal capacitance follows the spec:
%       C = (panel areal mass * A + equipment mass) * c_p

in = {A, coat.alpha, coat.eps_c, coat.eps_h, coat.Tt_C, coat.w, ...
      loadCase.Qint_sun, loadCase.Qint_ecl, p.sc.equip_mass};
n = cellfun(@numel, in);
N = max(n);
assert(all(n == 1 | n == N), 'starc:config:size', ...
    'Inputs must be scalars or share one length N.');
col = @(x) x(:) .* ones(N, 1);

cfg.N     = N;
cfg.A     = col(A);
cfg.C     = (p.sc.panel_areal_mass .* cfg.A + col(p.sc.equip_mass)) .* p.sc.cp;
cfg.alpha = col(coat.alpha);
cfg.eps_c = col(coat.eps_c);
cfg.eps_h = col(coat.eps_h);
cfg.Tt    = col(coat.Tt_C) + p.T0K;       % transition temperature [K]
cfg.w     = col(coat.w);                  % transition width [K]
cfg.Qsun  = col(loadCase.Qint_sun);       % internal dissipation, sunlit [W]
cfg.Qecl  = col(loadCase.Qint_ecl);       % internal dissipation, eclipse [W]
cfg.Tmin  = col(p.sc.T_min_C) + p.T0K;    % heater set point [K]

assert(all(cfg.A > 0) && all(cfg.w > 0), 'starc:config:value', ...
    'Area and transition width must be positive.');
end
