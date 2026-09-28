function e = starc_economics(E_ecl_J, E_sun_J, p)
%STARC_ECONOMICS  Battery + array mass and dollar value of heater energy.
%
%   e = starc_economics(E_ecl_J, E_sun_J, p)
%
%   E_ecl_J : heater energy drawn during eclipse, per orbit [J]
%   E_sun_J : heater energy drawn in sunlight, per orbit [J]
%   (arrays allowed; pass a SAVING to get the value of that saving)
%
%   Follows white paper Step 2 / Table 3:
%     battery  = eclipse energy / usable specific energy
%     array    = power needed over the sunlit period to (a) recharge the
%                eclipse energy through the recharge efficiency and
%                (b) run the sunlit heater directly, / array specific power
%   e.g. 1 W through a 35 min eclipse -> 3.9 g battery + 11.4 g array.
%
%   A "heater watt" in the paper is 1 W held through the whole eclipse.
%   W_equiv converts any heater energy profile into that unit by mass,
%   and the paper's $150-$800 per watt is applied to W_equiv.

ec = p.econ;
t_sun_h = p.orbit.sun_min / 60;
t_ecl_h = (p.orbit.period_min - p.orbit.sun_min) / 60;

E_ecl_Wh = E_ecl_J / 3600;
E_sun_Wh = E_sun_J / 3600;

e.battery_kg = E_ecl_Wh ./ ec.battery_Wh_per_kg;
e.array_W    = (E_ecl_Wh ./ ec.recharge_eff + E_sun_Wh) ./ t_sun_h;
e.array_kg   = e.array_W ./ ec.array_W_per_kg;
e.mass_kg    = e.battery_kg + e.array_kg;

e.battery_g = 1e3 * e.battery_kg;
e.array_g   = 1e3 * e.array_kg;
e.mass_g    = 1e3 * e.mass_kg;

% Mass for one eclipse heater watt (white paper: ~15.3 g)
e.kg_per_eclipse_W = t_ecl_h / ec.battery_Wh_per_kg + ...
                     (t_ecl_h / ec.recharge_eff) / t_sun_h / ec.array_W_per_kg;
e.W_equiv = e.mass_kg ./ e.kg_per_eclipse_W;

e.usd_launch = e.mass_kg .* ec.launch_usd_per_kg;   % launch cost of the mass only
e.usd_low    = e.W_equiv .* ec.usd_per_W(1);        % constellation value
e.usd_high   = e.W_equiv .* ec.usd_per_W(end);      % science-mission value
end
