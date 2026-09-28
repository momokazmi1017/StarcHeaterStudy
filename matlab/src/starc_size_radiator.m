function [A, info] = starc_size_radiator(p, coat, loadCase, dt)
%STARC_SIZE_RADIATOR  Radiator area whose hot-case orbit peaks at T_max.
%
%   [A, info] = starc_size_radiator(p, coat, loadCase, dt)
%
%   For each case, finds A such that the periodic steady-state peak
%   temperature of the hot load case equals p.sc.T_max_C. Peak temperature
%   falls monotonically with area (more emitting surface, more thermal
%   mass), so a bracketed Illinois (modified regula falsi) iteration on
%   log(A) is used, all cases in parallel.
%
%   coat and loadCase may hold column vectors (see STARC_CONFIG); A is N-by-1.

Tmax = p.sc.T_max_C + p.T0K;
cfg0 = starc_config(p, coat, loadCase, 1);
N = cfg0.N;
f = @(x) peak_error(exp(x), p, coat, loadCase, dt, Tmax);

a = log(p.num.A_bracket_m2(1)) * ones(N, 1);   fa = f(a);
b = log(p.num.A_bracket_m2(2)) * ones(N, 1);   fb = f(b);
assert(all(fa > 0) && all(fb < 0), 'starc:size:bracket', ...
    ['T_max cannot be met inside p.num.A_bracket_m2 for every case ' ...
     '(too hot even at the largest area, or too cold at the smallest).']);

x = b;  fx = fb;
for it = 1:p.num.size_maxit
    x  = (a .* fb - b .* fa) ./ (fb - fa);
    fx = f(x);
    if all(abs(fx) < p.num.size_tol_K), break; end
    sameB = sign(fx) == sign(fb);
    fa(sameB)  = fa(sameB) / 2;                  % Illinois down-weighting
    a(~sameB)  = b(~sameB);
    fa(~sameB) = fb(~sameB);
    b = x;  fb = fx;
end

A = exp(x);
info.iterations = it;
info.peak_error_K = fx;
info.converged = abs(fx) < p.num.size_tol_K;
end

function d = peak_error(A, p, coat, loadCase, dt, Tmax)
sim = starc_simulate(starc_config(p, coat, loadCase, A), p, dt, false);
d = sim.Tpeak_K - Tmax;
end
