function q = starc_radiated_flux(eps, T_C, p)
%STARC_RADIATED_FLUX  Grey-body emission to deep space, eps*sigma*T^4 [W/m^2].
%
%   q = starc_radiated_flux(eps, T_C, p)   T_C in degC, p from STARC_DEFAULT_PARAMS.

q = eps .* p.sigma .* (T_C + p.T0K).^4;
end
