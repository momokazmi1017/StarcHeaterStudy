function e = starc_emissivity(T, eps_c, eps_h, Tt, w)
%STARC_EMISSIVITY  Logistic emissivity model eps(T).
%
%   e = starc_emissivity(T, eps_c, eps_h, Tt, w)
%
%   eps(T) = eps_c + (eps_h - eps_c) / (1 + exp(-(T - Tt)/w))
%
%   T and Tt must share units (K or degC); w is in K. All inputs may be
%   arrays that broadcast. With eps_c == eps_h this is a fixed coating.
%
%   This is an idealised curve for a first-order study; the real STARC
%   response (shape, hysteresis, angle/wavelength dependence) is not
%   public and is not modelled.

e = eps_c + (eps_h - eps_c) ./ (1 + exp(-(T - Tt) ./ w));
end
