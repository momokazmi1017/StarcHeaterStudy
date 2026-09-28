function A = starc_datacenter_area(P_compute_W, overhead, T_K, eps, q_env, nFaces, sigma)
%STARC_DATACENTER_AREA  First-order radiator panel area for an orbital data center.
%
%   A = starc_datacenter_area(P_compute_W, overhead, T_K, eps, q_env, nFaces, sigma)
%
%   Heat to reject = P_compute * (1 + overhead). Each panel face radiates
%   eps*sigma*T^4 and absorbs q_env [W/m^2] from the environment, so
%
%       A = P (1 + overhead) / (nFaces * (eps sigma T^4 - q_env))
%
%   White paper example: 1 MW, 10 %, 330 K, eps 0.90, 150 W/m^2, two
%   faces -> about 1,200 m^2 of panel.

A = P_compute_W .* (1 + overhead) ./ (nFaces .* (eps .* sigma .* T_K.^4 - q_env));
end
