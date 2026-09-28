%RUN_STARC_STUDY  First-order heater-power study for an adaptive-emissivity radiator.
%
%   Run from anywhere:   >> run_starc_study
%
%   Compares a STARC-like adaptive-emissivity coating with fixed silver FEP
%   and glass OSR on a one-node LEO radiator:
%     1. sizes each coating's radiator so the hot case peaks at T_max,
%     2. runs the cold case on that area with a thermostatic heater,
%     3. converts the heater energy saved to battery + array mass and $,
%     4. sweeps the transition temperature T_t,
%     5. runs a 500-case Monte Carlo over uncertain inputs.
%   Plots go to results/*.png and key numbers to results/summary.txt.
%
%   All inputs live in src/starc_default_params.m.
%
%   FIRST-ORDER PUBLIC-DATA ESTIMATE. NOT VALIDATED AGAINST STARC TEST DATA.
%
%   Requires MATLAB R2021a or later, no toolboxes.

here   = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
outDir = fullfile(here, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

p     = starc_default_params();
names = {'OSR', 'FEP', 'STARC'};           % fixed order = fixed plot colours
iS    = find(strcmp(names, 'STARC'));
iB    = find(strcmp(names, p.baseline));
iFix  = find(~strcmp(names, 'STARC'));
sty   = plot_style(names);
t0    = tic;

fprintf(['\nSTARC heater-power study\n' ...
         '  First-order public-data estimate - NOT validated against STARC data.\n\n']);

%% 1. Nominal fair comparison (headline results, dt = p.num.dt)
fprintf('[1/5] Sizing radiators and running hot/cold cases (dt = %g s) ...\n', p.num.dt);
R   = starc_compare(p, names, p.num.dt);
sav = coating_savings(R, iS, p);

% Check that the coarser step used for the sweep and Monte Carlo agrees
Rs = starc_compare(p, names, p.num.dt_study);
Escale = max(max(R.cold.E_heater_J), 1);
dtCheck.A_rel = max(abs(Rs.A - R.A) ./ R.A);
dtCheck.E_rel = max(abs(Rs.cold.E_heater_J - R.cold.E_heater_J)) / Escale;

%% 2. Plots 1-3: temperature, heater power, eps(T)
fprintf('[2/5] Plotting orbit histories and eps(T) ...\n');
plot_temperature(R, p, sty, fullfile(outDir, 'fig1_temperature_orbit.png'));
plot_heater(R, p, sty, fullfile(outDir, 'fig2_heater_power_cold.png'));
plot_emissivity(p, sty, names, fullfile(outDir, 'fig3_starc_emissivity.png'));

%% 3. Transition-temperature sweep (plot 4)
fprintf('[3/5] Sweeping transition temperature T_t (dt = %g s) ...\n', p.num.dt_study);
sw.Tt_C = p.study.Tt_sweep_C(:);
cS = p.coat.STARC;  cS.Tt_C = sw.Tt_C;
sw.A    = starc_size_radiator(p, cS, p.hot, p.num.dt_study);
sw.cold = starc_simulate(starc_config(p, cS, p.cold, sw.A), p, p.num.dt_study, false);
% Heater energy saved vs each fixed coating [Wh/orbit], columns follow iFix
sw.saved_Wh = (Rs.cold.E_heater_J(iFix).' - sw.cold.E_heater_J) / 3600;
% Best switch point. "Good" T_t values save within 0.5 % (of the baseline
% heater energy) of the maximum AND need a radiator within 0.1 % of the
% smallest such area. Often this is a window, not a point; pick its centre
% so the choice has the most margin against T_t uncertainty.
sb   = sw.saved_Wh(:, iFix == iB);
tolW = 0.005 * Rs.cold.E_heater_J(iB) / 3600;
cand = find(sb >= max(sb) - tolW);
good = cand(sw.A(cand) <= 1.001 * min(sw.A(cand)));
sw.window_C = [min(sw.Tt_C(good)), max(sw.Tt_C(good))];
[~, m] = min(abs(sw.Tt_C(good) - mean(sw.window_C)));
sw.iBest = good(m);
plot_sweep(sw, Rs, p, sty, iFix, iB, iS, fullfile(outDir, 'fig4_transition_temperature_sweep.png'));

%% 4. Monte Carlo (plot 5)
fprintf('[4/5] Monte Carlo, %d runs vs %s (dt = %g s) ...\n', p.mc.N, p.baseline, p.num.dt_study);
mc = run_monte_carlo(p);
plot_mc(mc, p, sty, iS, fullfile(outDir, 'fig5_monte_carlo_heater_watts_saved.png'));

%% 5. Summary
fprintf('[5/5] Writing summary ...\n\n');
txt = build_summary(p, R, sav, names, iS, iB, iFix, sw, mc, dtCheck, toc(t0));
fid = fopen(fullfile(outDir, 'summary.txt'), 'w');
fprintf(fid, '%s', txt);
fclose(fid);
fprintf('%s', txt);
fprintf('\nResults written to %s\n', outDir);


%% ========================================================================
%  Local functions
%  ========================================================================

function sav = coating_savings(R, iS, p)
% STARC savings relative to every coating in R (row iS is STARC itself)
ft2 = p.econ.ft2_per_m2;
for j = 1:numel(R.A)
    dEecl = R.cold.E_heater_ecl_J(j) - R.cold.E_heater_ecl_J(iS);
    dEsun = R.cold.E_heater_sun_J(j) - R.cold.E_heater_sun_J(iS);
    s.E_saved_Wh     = (dEecl + dEsun) / 3600;
    s.Qavg_saved_W   = R.cold.Qh_avg_W(j)  - R.cold.Qh_avg_W(iS);
    s.Qpeak_saved_W  = R.cold.Qh_peak_W(j) - R.cold.Qh_peak_W(iS);
    s.econ           = starc_economics(dEecl, dEsun, p);
    s.area_saved_pct = 100 * (R.A(j) - R.A(iS)) / R.A(j);
    s.usd_per_ft2    = [s.econ.usd_low, s.econ.usd_high] / (R.A(iS) * ft2);
    s.Wequiv_per_ft2 = s.econ.W_equiv / (R.A(iS) * ft2);
    sav(j) = s; %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function mc = run_monte_carlo(p)
% Paired Monte Carlo: each draw is applied to STARC and to the baseline
% coating (same loads and alpha change), both re-sized, then compared.
rng(p.mc.seed);
N = p.mc.N;
U = @(r) r(1) + (r(2) - r(1)) * rand(N, 1);
d.w      = U(p.mc.w_K);
d.eps_c  = U(p.mc.eps_c);
d.eps_h  = U(p.mc.eps_h);
d.kHot   = U(p.mc.Qint_scale);
d.kCold  = U(p.mc.Qint_scale);
d.dalpha = U(p.mc.dalpha);

hot.Qint_sun  = p.hot.Qint_sun  * d.kHot;   hot.Qint_ecl  = p.hot.Qint_ecl  * d.kHot;
cold.Qint_sun = p.cold.Qint_sun * d.kCold;  cold.Qint_ecl = p.cold.Qint_ecl * d.kCold;

cS = p.coat.STARC;
cS.w = d.w;  cS.eps_c = d.eps_c;  cS.eps_h = d.eps_h;  cS.alpha = cS.alpha + d.dalpha;
cB = p.coat.(p.baseline);
cB.alpha = cB.alpha + d.dalpha;

dt = p.num.dt_study;
mc.A_S = starc_size_radiator(p, cS, hot, dt);
mc.A_B = starc_size_radiator(p, cB, hot, dt);
colS = starc_simulate(starc_config(p, cS, cold, mc.A_S), p, dt, false);
colB = starc_simulate(starc_config(p, cB, cold, mc.A_B), p, dt, false);

mc.econ = starc_economics(colB.E_heater_ecl_J - colS.E_heater_ecl_J, ...
                          colB.E_heater_sun_J - colS.E_heater_sun_J, p);
mc.W_saved        = mc.econ.W_equiv;
mc.area_saved_pct = 100 * (mc.A_B - mc.A_S) ./ mc.A_B;
mc.pW   = pctl(mc.W_saved, [10 50 90]);
mc.pA   = pctl(mc.area_saved_pct, [10 50 90]);
mc.pUSD = pctl(mc.econ.usd_low, [10 50 90]);
mc.pUSDhi = pctl(mc.econ.usd_high, [10 50 90]);
mc.draws = d;
mc.converged = all(colS.converged) && all(colB.converged);
end

% -------------------------------------------------------------------------
function q = pctl(x, pc)
% Percentiles with linear interpolation (no Statistics Toolbox needed)
x   = sort(x(:));
pos = 1 + (numel(x) - 1) * pc(:).' / 100;
lo  = floor(pos);  hi = ceil(pos);
xl  = x(lo);  xh = x(hi);
q   = xl(:).' + (pos - lo) .* (xh(:).' - xl(:).');
end

% -------------------------------------------------------------------------
function s = plot_style(names)
% Fixed colour per coating (validated categorical slots 1-3) plus a
% distinct line style so identity never relies on colour alone.
pal = struct('OSR',   [42 120 214] / 255, ...   % blue
             'FEP',   [235 104 52] / 255, ...   % orange
             'STARC', [27 175 122] / 255);      % aqua
lsty = struct('OSR', '--', 'FEP', '-.', 'STARC', '-');
for j = 1:numel(names)
    s.color(j, :) = pal.(names{j});
    s.line{j}     = lsty.(names{j});
end
s.ink   = [0.32 0.32 0.31];     % secondary text / reference lines
s.shade = [0.55 0.55 0.55];
s.deg   = [char(176) 'C'];
end

function f = new_figure(p, sz)
vis = 'off';
if p.out.show_figures, vis = 'on'; end
f = figure('Color', 'w', 'Visible', vis, 'Position', [100 100 sz]);
end

function style_axes(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10, ...
    'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3], ...
    'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1);
grid(ax, 'on');
end

function shade_eclipse(ax, p, sty)
% Grey band over the eclipse part of the orbit
yl = ylim(ax);
x1 = p.orbit.sun_min;  x2 = p.orbit.period_min;
h = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], sty.shade, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(h, 'bottom');
text(ax, (x1 + x2) / 2, yl(2), 'eclipse', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'Color', sty.ink, 'FontSize', 9);
ylim(ax, yl);
xlim(ax, [0 p.orbit.period_min]);
end

% -------------------------------------------------------------------------
function plot_temperature(R, p, sty, file)
f  = new_figure(p, [900 680]);
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
cases  = {R.hot, R.cold};
ttl    = {sprintf('Hot case: Q_{int} = %g W, radiator sized so peak = %g %s', ...
                  p.hot.Qint_sun, p.sc.T_max_C, sty.deg), ...
          sprintf('Cold case: Q_{int} = %g W on the same area, heater holds T \\geq %g %s', ...
                  p.cold.Qint_sun, p.sc.T_min_C, sty.deg)};
limVal = [p.sc.T_max_C, p.sc.T_min_C];
limTxt = {'T_{max}', 'T_{min}'};
for c = 1:2
    ax = nexttile(tl);  hold(ax, 'on');  style_axes(ax);
    sim = cases{c};
    for j = 1:numel(R.names)
        plot(ax, sim.t_min, sim.T_C(j, :), 'Color', sty.color(j, :), ...
            'LineStyle', sty.line{j}, 'LineWidth', 2, ...
            'DisplayName', sprintf('%s  (A = %.2f m^2)', R.names{j}, R.A(j)));
    end
    yline(ax, limVal(c), ':', limTxt{c}, 'Color', sty.ink, 'LineWidth', 1, ...
        'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'left');
    yl = ylim(ax);  ylim(ax, [min(yl(1), limVal(c) - 3), max(yl(2), limVal(c) + 3)]);
    shade_eclipse(ax, p, sty);
    title(ax, ttl{c}, 'FontWeight', 'normal');
    ylabel(ax, ['Radiator temperature (' sty.deg ')']);
    legend(ax, 'Location', 'eastoutside', 'Box', 'off');
end
xlabel(ax, 'Time in orbit (min)');
exportgraphics(f, file, 'Resolution', 150);
end

% -------------------------------------------------------------------------
function plot_heater(R, p, sty, file)
f  = new_figure(p, [900 420]);
ax = axes(f);  hold(ax, 'on');  style_axes(ax);
for j = 1:numel(R.names)
    plot(ax, R.cold.tq_min, R.cold.Qh_W(j, :), 'Color', sty.color(j, :), ...
        'LineStyle', sty.line{j}, 'LineWidth', 2, ...
        'DisplayName', sprintf('%s: %.1f Wh/orbit, peak %.0f W', R.names{j}, ...
                               R.cold.E_heater_J(j) / 3600, R.cold.Qh_peak_W(j)));
end
yl = ylim(ax);  ylim(ax, [0 max(yl(2), 1)]);
shade_eclipse(ax, p, sty);
title(ax, sprintf('Heater power, cold case (Q_{int} = %g W, T_{min} = %g %s)', ...
      p.cold.Qint_sun, p.sc.T_min_C, sty.deg), 'FontWeight', 'normal');
xlabel(ax, 'Time in orbit (min)');
ylabel(ax, 'Heater power (W)');
legend(ax, 'Location', 'eastoutside', 'Box', 'off');
exportgraphics(f, file, 'Resolution', 150);
end

% -------------------------------------------------------------------------
function plot_emissivity(p, sty, names, file)
f  = new_figure(p, [800 460]);
ax = axes(f);  hold(ax, 'on');  style_axes(ax);
T  = -60:0.25:80;
for j = 1:numel(names)
    c = p.coat.(names{j});
    e = starc_emissivity(T, c.eps_c, c.eps_h, c.Tt_C, c.w);
    lbl = c.name;
    if c.eps_c ~= c.eps_h
        lbl = sprintf('%s  (T_t = %g %s, w = %g K)', c.name, c.Tt_C, sty.deg, c.w);
    end
    plot(ax, T, e, 'Color', sty.color(j, :), 'LineStyle', sty.line{j}, ...
        'LineWidth', 2, 'DisplayName', lbl);
end
xline(ax, p.sc.T_min_C, ':', 'T_{min}', 'Color', sty.ink, 'HandleVisibility', 'off');
xline(ax, p.sc.T_max_C, ':', 'T_{max}', 'Color', sty.ink, 'HandleVisibility', 'off');
ylim(ax, [0 1]);
title(ax, 'Infrared emissivity vs temperature (idealised logistic STARC curve)', ...
      'FontWeight', 'normal');
xlabel(ax, ['Surface temperature (' sty.deg ')']);
ylabel(ax, 'Emissivity \epsilon (-)');
legend(ax, 'Location', 'southeast', 'Box', 'off');
exportgraphics(f, file, 'Resolution', 150);
end

% -------------------------------------------------------------------------
function plot_sweep(sw, Rs, p, sty, iFix, iB, iS, file)
f  = new_figure(p, [860 680]);
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
b  = sw.iBest;

ax1 = nexttile(tl);  hold(ax1, 'on');  style_axes(ax1);
for k = 1:numel(iFix)
    j = iFix(k);
    plot(ax1, sw.Tt_C, sw.saved_Wh(:, k), 'Color', sty.color(j, :), ...
        'LineStyle', sty.line{j}, 'LineWidth', 2, ...
        'DisplayName', ['vs ' Rs.names{j}]);
end
kB = find(iFix == iB);
plot(ax1, sw.Tt_C(b), sw.saved_Wh(b, kB), 'o', 'MarkerSize', 9, 'LineWidth', 2, ...
    'MarkerEdgeColor', 'w', 'MarkerFaceColor', sty.color(iS, :), ...
    'DisplayName', sprintf('best T_t = %g %s', sw.Tt_C(b), sty.deg));
yline(ax1, 0, '-', 'Color', sty.ink, 'HandleVisibility', 'off');
yl = ylim(ax1);  win = sw.window_C;
h = patch(ax1, [win(1) win(2) win(2) win(1)], [yl(1) yl(1) yl(2) yl(2)], sty.color(iS, :), ...
    'FaceAlpha', 0.10, 'EdgeColor', 'none', 'DisplayName', 'good T_t window');
uistack(h, 'bottom');
ylim(ax1, yl);
title(ax1, 'Cold-case heater energy saved by STARC (each coating sized for its own hot case)', ...
      'FontWeight', 'normal');
ylabel(ax1, 'Heater energy saved (Wh/orbit)');
legend(ax1, 'Location', 'best', 'Box', 'off');

ax2 = nexttile(tl);  hold(ax2, 'on');  style_axes(ax2);
plot(ax2, sw.Tt_C, sw.A, 'Color', sty.color(iS, :), 'LineStyle', sty.line{iS}, ...
    'LineWidth', 2, 'DisplayName', 'STARC');
for k = 1:numel(iFix)
    j = iFix(k);
    yline(ax2, Rs.A(j), sty.line{j}, Rs.names{j}, 'Color', sty.color(j, :), ...
        'LineWidth', 1.5, 'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'left');
end
plot(ax2, sw.Tt_C(b), sw.A(b), 'o', 'MarkerSize', 9, 'LineWidth', 2, ...
    'MarkerEdgeColor', 'w', 'MarkerFaceColor', sty.color(iS, :), 'HandleVisibility', 'off');
xline(ax2, sw.Tt_C(b), ':', 'Color', sty.ink, 'HandleVisibility', 'off');
xline(ax1, sw.Tt_C(b), ':', 'Color', sty.ink, 'HandleVisibility', 'off');
title(ax2, 'Radiator area needed to hold the hot-case peak at T_{max}', 'FontWeight', 'normal');
ylabel(ax2, 'Radiator area (m^2)');
xlabel(ax2, ['STARC transition temperature T_t (' sty.deg ')']);
linkaxes([ax1 ax2], 'x');
xlim(ax2, [min(sw.Tt_C) max(sw.Tt_C)]);
exportgraphics(f, file, 'Resolution', 150);
end

% -------------------------------------------------------------------------
function plot_mc(mc, p, sty, iS, file)
f  = new_figure(p, [820 460]);
ax = axes(f);  hold(ax, 'on');  style_axes(ax);
histogram(ax, mc.W_saved, 30, 'FaceColor', sty.color(iS, :), 'FaceAlpha', 0.9, ...
    'EdgeColor', 'w', 'LineWidth', 1);
lbl = {'P10', 'P50', 'P90'};
for k = 1:3
    xline(ax, mc.pW(k), '--', sprintf('%s = %.1f W', lbl{k}, mc.pW(k)), ...
        'Color', sty.ink, 'LineWidth', 1.2, 'LabelVerticalAlignment', 'top');
end
title(ax, sprintf('Monte Carlo (%d runs): heater watts saved by STARC vs %s', ...
      p.mc.N, p.coat.(p.baseline).name), 'FontWeight', 'normal');
xlabel(ax, 'Heater power saved (equivalent eclipse watts, W)');
ylabel(ax, 'Runs');
exportgraphics(f, file, 'Resolution', 150);
end

% -------------------------------------------------------------------------
function txt = build_summary(p, R, sav, names, iS, iB, iFix, sw, mc, dtCheck, runtime)
ft2 = p.econ.ft2_per_m2;
d   = 'C';
L = {};
L{end+1} = '=====================================================================';
L{end+1} = ' STARC adaptive-emissivity radiator: first-order heater-power study';
L{end+1} = '=====================================================================';
L{end+1} = sprintf(' Generated %s', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
L{end+1} = '';
L{end+1} = ' *** FIRST-ORDER PUBLIC-DATA ESTIMATE. NOT VALIDATED AGAINST STARC ***';
L{end+1} = ' *** TEST DATA. One-node model; see README.md for all limitations. ***';
L{end+1} = '';
L{end+1} = ' SCENARIO';
L{end+1} = sprintf('   LEO orbit %g min (%g min sun, %g min eclipse); F_sun %.2f, F_earth %.2f', ...
    p.orbit.period_min, p.orbit.sun_min, p.orbit.period_min - p.orbit.sun_min, ...
    p.env.F_sun, p.env.F_earth);
L{end+1} = sprintf('   Hot case  Q_int %g W (sun) / %g W (eclipse); radiator sized so peak = %g %s', ...
    p.hot.Qint_sun, p.hot.Qint_ecl, p.sc.T_max_C, d);
L{end+1} = sprintf('   Cold case Q_int %g W (sun) / %g W (eclipse); heater holds T >= %g %s', ...
    p.cold.Qint_sun, p.cold.Qint_ecl, p.sc.T_min_C, d);
c = p.coat.STARC;
L{end+1} = sprintf('   STARC eps %.2f -> %.2f, T_t = %g %s, w = %g K, alpha %.2f', ...
    c.eps_c, c.eps_h, c.Tt_C, d, c.w, c.alpha);
L{end+1} = '';
L{end+1} = sprintf(' 1. RADIATOR AREA AND COLD-CASE HEATER DEMAND (dt = %g s)', p.num.dt);
L{end+1} = sprintf('   %-11s %8s %8s %12s %11s %10s %10s  %s', 'Coating', 'Area m2', ...
    'Area ft2', 'Heater Wh/o', 'eclipse Wh', 'Peak W', 'Avg W', 'Cold T range');
for j = 1:numel(names)
    L{end+1} = sprintf('   %-11s %8.3f %8.2f %12.2f %11.2f %10.1f %10.1f  %.1f to %.1f %s', ...
        R.names{j}, R.A(j), R.A(j) * ft2, R.cold.E_heater_J(j) / 3600, ...
        R.cold.E_heater_ecl_J(j) / 3600, R.cold.Qh_peak_W(j), R.cold.Qh_avg_W(j), ...
        R.cold.Tlow_C(j), R.cold.Tpeak_C(j), d); %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = ' 2. WHAT STARC SAVES';
for j = iFix
    s = sav(j);  e = s.econ;
    L{end+1} = sprintf('   vs %s%s', R.names{j}, ternary(j == iB, '  (baseline)', '')); %#ok<AGROW>
    L{end+1} = sprintf('     Heater energy saved        %8.1f Wh per orbit', s.E_saved_Wh); %#ok<AGROW>
    L{end+1} = sprintf('     Heater power saved         %8.1f W orbit-average, %.1f W peak', ...
        s.Qavg_saved_W, s.Qpeak_saved_W); %#ok<AGROW>
    L{end+1} = sprintf('     = equivalent eclipse watts %8.1f W  (paper Step 2 unit)', e.W_equiv); %#ok<AGROW>
    L{end+1} = sprintf('     Battery mass saved         %8.0f g', e.battery_g); %#ok<AGROW>
    L{end+1} = sprintf('     Solar-array mass saved     %8.0f g', e.array_g); %#ok<AGROW>
    L{end+1} = sprintf('     Total mass saved           %8.0f g  (launch only: $%s at $%s/kg)', ...
        e.mass_g, commas(e.usd_launch), commas(p.econ.launch_usd_per_kg)); %#ok<AGROW>
    L{end+1} = sprintf('     Value at $%g-$%g per W    $%s - $%s per spacecraft', ...
        p.econ.usd_per_W(1), p.econ.usd_per_W(end), commas(e.usd_low), commas(e.usd_high)); %#ok<AGROW>
    L{end+1} = sprintf('     Per sq ft of STARC radiator $%s - $%s  (%.1f W/sq ft; paper credits 3 W/sq ft)', ...
        commas(s.usd_per_ft2(1)), commas(s.usd_per_ft2(2)), s.Wequiv_per_ft2); %#ok<AGROW>
    L{end+1} = sprintf('     Radiator area saved        %8.1f %%', s.area_saved_pct); %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = ' 3. TRANSITION-TEMPERATURE SWEEP';
b = sw.iBest;  kB = find(iFix == iB);
L{end+1} = sprintf('   Best T_t for this mission: %g %s  (saves %.1f Wh/orbit vs %s, area %.3f m2).', ...
    sw.Tt_C(b), d, sw.saved_Wh(b, kB), R.names{iB}, sw.A(b));
L{end+1} = sprintf(['   Good window: T_t = %g to %g %s (heater saving within 0.5%% of max and\n' ...
    '   radiator within 0.1%% of the smallest); best = centre of that window.'], ...
    sw.window_C(1), sw.window_C(2), d);
L{end+1} = '   Too low a T_t keeps eps high when cold (heater still needed); too high';
L{end+1} = '   a T_t leaves eps low when hot, so the radiator must grow.';
L{end+1} = '';
L{end+1} = sprintf(' 4. MONTE CARLO (%d paired runs vs %s)', p.mc.N, R.names{iB});
L{end+1} = sprintf('   Uncertain: w %g-%g K, eps_c %.2f-%.2f, eps_h %.2f-%.2f, Q_int x%.1f-%.1f,', ...
    p.mc.w_K, p.mc.eps_c, p.mc.eps_h, p.mc.Qint_scale);
L{end+1} = sprintf('   alpha change %+.2f to %+.2f (uniform distributions).', p.mc.dalpha);
L{end+1} = sprintf('   Heater watts saved (eclipse-W)   P10 %.1f   P50 %.1f   P90 %.1f W', mc.pW);
L{end+1} = sprintf('   Value at $%g/W                  P10 $%s  P50 $%s  P90 $%s', ...
    p.econ.usd_per_W(1), commas(mc.pUSD(1)), commas(mc.pUSD(2)), commas(mc.pUSD(3)));
L{end+1} = sprintf('   Value at $%g/W                  P10 $%s  P50 $%s  P90 $%s', ...
    p.econ.usd_per_W(end), commas(mc.pUSDhi(1)), commas(mc.pUSDhi(2)), commas(mc.pUSDhi(3)));
L{end+1} = sprintf('   Radiator area saved              P10 %.1f%%  P50 %.1f%%  P90 %.1f%%', mc.pA);
L{end+1} = '';
L{end+1} = ' 5. NUMERICAL CHECKS';
L{end+1} = sprintf('   Periodic energy closure, last orbit: max |in-out|/out = %.1e (hot), %.1e (cold)', ...
    max(abs(R.hot.closure)), max(abs(R.cold.closure)));
L{end+1} = sprintf('   dt = %g s vs %g s: area differs by %.3f%%, heater energy by %.3f%% of max', ...
    p.num.dt_study, p.num.dt, 100 * dtCheck.A_rel, 100 * dtCheck.E_rel);
L{end+1} = sprintf('   Periodic solves converged: %s;  sizing converged: %s;  Monte Carlo: %s', ...
    yesno(all(R.hot.converged) && all(R.cold.converged)), yesno(all(R.sizing.converged)), ...
    yesno(mc.converged));
L{end+1} = sprintf('   Run time %.0f s', runtime);
L{end+1} = '';
L{end+1} = ' Heater watts are expressed in the white paper''s unit: the battery + array';
L{end+1} = ' mass of 1 W held through the 35 min eclipse (~15.3 g). $ values apply the';
L{end+1} = ' paper''s $150 (constellation) to $800 (science) per watt; they are gross';
L{end+1} = ' savings, before the price difference of STARC vs the incumbent coating.';
txt = [strjoin(L, newline) newline];
end

function s = commas(x)
% Round to a readable precision and add thousands separators
if abs(x) >= 100, x = round(x, -1); else, x = round(x); end
s = sprintf('%.0f', abs(x));
s = regexprep(s, '(\d)(?=(\d{3})+$)', '$1,');
if x < 0, s = ['-' s]; end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function s = yesno(c)
s = ternary(c, 'yes', 'NO');
end
