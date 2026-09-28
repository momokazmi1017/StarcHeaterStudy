function R = starc_compare(p, names, dt)
%STARC_COMPARE  Fair comparison of coatings: size on the hot case, run the cold case.
%
%   R = starc_compare(p, names, dt)      e.g. names = {'FEP','OSR','STARC'}
%
%   1. Each coating gets its own radiator area so the hot case (sunlit,
%      high Q_int) peaks exactly at T_max.
%   2. The cold case (low Q_int) is run on that same area with the
%      thermostatic heater holding T >= T_min.
%
%   R.names, R.keys : coating display names / keys in p.coat
%   R.A             : radiator area [m^2]
%   R.hot, R.cold   : STARC_SIMULATE outputs (with time histories)

if nargin < 3, dt = p.num.dt; end
coat = starc_stack_coatings(p, names);

[R.A, R.sizing] = starc_size_radiator(p, coat, p.hot, dt);
R.hot  = starc_simulate(starc_config(p, coat, p.hot,  R.A), p, dt, true);
R.cold = starc_simulate(starc_config(p, coat, p.cold, R.A), p, dt, true);
R.names = coat.name;
R.keys  = coat.key;
R.dt    = dt;
end
