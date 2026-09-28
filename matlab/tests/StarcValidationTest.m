classdef StarcValidationTest < matlab.unittest.TestCase
%STARCVALIDATIONTEST  Checks the model against the white paper and itself.
%
%   Run all tests from the project folder:
%       results = runtests('tests');  table(results)
%
%   Hand-calc checks reproduce numbers printed in Aetara White Paper 01;
%   the remaining tests check the solver (limit case, energy closure,
%   time-step convergence). None of them validate the model against
%   STARC hardware data, which is not public.

    properties
        p   % default parameters
    end

    methods (TestClassSetup)
        function setup(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
            tc.p = starc_default_params();
        end
    end

    methods (Test)
        function step1RadiatedFlux(tc)
            % White paper Step 1: at -10 C, eps 0.80 radiates ~217 W/m^2 and
            % eps 0.20 ~54 W/m^2.
            tc.verifyEqual(starc_radiated_flux(0.80, -10, tc.p), 217, 'AbsTol', 2);
            tc.verifyEqual(starc_radiated_flux(0.20, -10, tc.p), 54,  'AbsTol', 1);
        end

        function heaterWattCost(tc)
            % White paper Step 2: 1 W through a 35 min eclipse ->
            % ~3.9 g battery + ~11.4 g array = ~15.3 g.
            e = starc_economics(1 * 35 * 60, 0, tc.p);
            tc.verifyEqual(e.battery_g, 3.9,  'AbsTol', 0.05);
            tc.verifyEqual(e.array_g,   11.4, 'AbsTol', 0.1);
            tc.verifyEqual(e.mass_g,    15.3, 'AbsTol', 0.1);
            tc.verifyEqual(e.W_equiv,   1,    'RelTol', 1e-12);
        end

        function orbitalDataCenterArea(tc)
            % White paper: 330 K, eps 0.90, two faces, 150 W/m^2 per face,
            % 1 MW compute + 10 % overhead -> ~1,200 m^2 of panel.
            A = starc_datacenter_area(1e6, 0.10, 330, 0.90, 150, 2, tc.p.sigma);
            tc.verifyEqual(A, 1200, 'RelTol', 0.02);
        end

        function logisticCurveShape(tc)
            % eps(T_t) is the midpoint; far from T_t it tends to eps_c / eps_h.
            c = tc.p.coat.STARC;
            f = @(T) starc_emissivity(T, c.eps_c, c.eps_h, c.Tt_C, c.w);
            tc.verifyEqual(f(c.Tt_C), (c.eps_c + c.eps_h) / 2, 'AbsTol', 1e-12);
            tc.verifyEqual(f(c.Tt_C - 20 * c.w), c.eps_c, 'AbsTol', 1e-6);
            tc.verifyEqual(f(c.Tt_C + 20 * c.w), c.eps_h, 'AbsTol', 1e-6);
            tc.verifyTrue(all(diff(f(-50:50)) > 0));
        end

        function equalEmissivitiesGiveFixedCoating(tc)
            % STARC with eps_c = eps_h = 0.88 (same alpha) must behave
            % exactly like glass OSR: same sized area and heater demand.
            p = tc.p;
            osr = p.coat.OSR;
            s = p.coat.STARC;  s.eps_c = osr.eps_c;  s.eps_h = osr.eps_h;  s.alpha = osr.alpha;
            coat = struct('alpha', [osr.alpha; s.alpha], 'eps_c', [osr.eps_c; s.eps_c], ...
                          'eps_h', [osr.eps_h; s.eps_h], 'Tt_C', [osr.Tt_C; s.Tt_C], ...
                          'w', [osr.w; s.w]);
            dt = p.num.dt_study;
            A = starc_size_radiator(p, coat, p.hot, dt);
            cold = starc_simulate(starc_config(p, coat, p.cold, A), p, dt, true);
            tc.verifyEqual(A(2), A(1), 'RelTol', 1e-9);
            tc.verifyEqual(cold.E_heater_J(2), cold.E_heater_J(1), 'RelTol', 1e-9);
            tc.verifyEqual(cold.T_C(2, :), cold.T_C(1, :), 'AbsTol', 1e-9);
            e = starc_emissivity(-50:80, s.eps_c, s.eps_h, s.Tt_C, s.w);
            tc.verifyEqual(e, 0.88 * ones(size(e)), 'AbsTol', 1e-15);
        end

        function periodicEnergyClosure(tc)
            % Over the reported (last) orbit, energy in = energy out within 1 %.
            R = starc_compare(tc.p, {'FEP', 'OSR', 'STARC'}, tc.p.num.dt);
            tc.verifyTrue(all(R.hot.converged) && all(R.cold.converged));
            tc.verifyLessThan(max(abs(R.hot.closure)),  0.01);
            tc.verifyLessThan(max(abs(R.cold.closure)), 0.01);
            % Independent check: recompute radiated energy from the recorded
            % temperature history with the trapezoid rule.
            p = tc.p;
            cfg = starc_config(p, starc_stack_coatings(p, R.keys), p.cold, R.A);
            TK = R.cold.T_K;
            e = starc_emissivity(TK, cfg.eps_c, cfg.eps_h, cfg.Tt, cfg.w);
            Eout = trapz(R.cold.t_min * 60, e .* p.sigma .* cfg.A .* TK.^4, 2);
            tc.verifyEqual(Eout, R.cold.E_out_J, 'RelTol', 1e-3);
            % Sized hot case really peaks at T_max
            tc.verifyEqual(R.hot.Tpeak_C, p.sc.T_max_C * ones(3, 1), 'AbsTol', 0.01);
        end

        function timestepConvergence(tc)
            % Halving dt (1 s -> 0.5 s) changes areas, heater energy and
            % STARC's heater saving by less than 0.5 %.
            names = {'FEP', 'OSR', 'STARC'};
            R1 = starc_compare(tc.p, names, 1);
            R2 = starc_compare(tc.p, names, 0.5);
            tc.verifyLessThan(max(abs(R2.A - R1.A) ./ R1.A), 0.005);
            Escale = max(max(R1.cold.E_heater_J), 1);
            tc.verifyLessThan(max(abs(R2.cold.E_heater_J - R1.cold.E_heater_J)) / Escale, 0.005);
            s1 = R1.cold.E_heater_J(1:2) - R1.cold.E_heater_J(3);
            s2 = R2.cold.E_heater_J(1:2) - R2.cold.E_heater_J(3);
            tc.verifyLessThan(max(abs(s2 - s1)) / max(max(abs(s1)), 1), 0.005);
        end
    end
end
