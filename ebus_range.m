function varargout = ebus_range(cmd, varargin)
%EBUS_RANGE Electric bus drive-cycle energy and range model.
%
%   ebus_range all        every cycle, one parameter set, one table (+ figures)
%   ebus_range sweeps     parameter sweeps: passengers, grade, aux, regen, speed
%   ebus_range tests      the validation suite
%   ebus_range profiles   export pack-current CSVs for a battery model
%   ebus_range help       this text
%
%   R = ebus_range('all') also returns the per-cycle result structs.
%
%   The whole project is this one file: the road-load model, the powertrain,
%   the five published drive cycles, the plots and the tests. Nothing to add
%   to the path and nothing to install -- just run it.
%
%   Figures land in ./figures and exported profiles in ./profiles, both
%   created relative to wherever you run it from.

if nargin < 1 || isempty(cmd); cmd = 'help'; end
varargout = {};

switch lower(cmd)
    case {'all', 'run_all'}
        if nargout > 0; varargout{1} = run_all(); else; run_all(); end
    case {'sweeps', 'run_sweeps'}
        run_sweeps();
    case {'tests', 'run_tests', 'test'}
        run_tests();
    case {'profiles', 'export_profiles'}
        export_profiles();
    case {'help', '-h', '--help'}
        disp(help(mfilename));
    otherwise
        error('ebus_range:cmd', ...
              'Unknown command "%s". Try: all, sweeps, tests, profiles, help.', cmd);
end
end


% ===========================================================================
% Model
% ===========================================================================

function p = params_bus(varargin)
%PARAMS_BUS Parameters for a 40 ft battery-electric transit bus.
%   P = PARAMS_BUS() returns the baseline (summer) parameter set.
%   P = PARAMS_BUS('n_pass', 60) overrides any field by name.
%
%   Representative figures from transit-electrification literature, not any
%   OEM's data. Every value has its basis written next to it, because a range
%   number is only as arguable as the parameters behind it.

% environment
p.g   = 9.80665;           % m/s^2
p.rho = 1.225;             % kg/m^3, air density at 15 C

% mass
p.m_curb  = 15000;         % kg, typical 40 ft BEB (~33,000 lb; GVWR ~45,000 lb)
p.n_pass  = 40;            % passengers: 38 seated + a few standing
p.m_pass  = 68;            % kg each, transit convention
p.eps_rot = 0.05;          % rotating inertia (wheels, gearbox, rotor) as a
                           % fraction of mass, added only to the accel term

% road load
p.Cd    = 0.60;            % blunt bus body
p.A     = 8.0;             % m^2 frontal area (2.6 m wide x 3.2 m tall)
p.Crr   = 0.008;           % radial truck tyres, dry pavement
p.grade = 0;               % rise/run, scalar or per-sample vector

% powertrain
p.eta_dt         = 0.90;   % battery -> wheel: inverter x motor x gearbox
p.eta_regen      = 0.85;   % wheel -> battery on the regen path
p.regen_enabled  = true;
p.P_regen_max    = 250e3;  % W, motor's rated regen capability
p.P_mot_max      = 250e3;  % W, rated traction power (feasibility check only)
p.v_regen_cutoff = 5/3.6;  % m/s; below this the friction brakes take over

% auxiliaries: HVAC, air compressor, steering, lighting. Drawn whether or not
% the bus is moving, which is why slow routes pay for them twice.
p.P_aux = 5e3;             % W, mild weather

% battery
p.capacity_kwh = 440;      % nameplate pack energy
p.soc_window   = 0.90;     % usable fraction (fleets do not run 100-0 %)
p.V_pack_nom   = 650;      % V, nominal (used only for the bms-sim export)

p.label = 'baseline';

% name/value overrides
assert(mod(numel(varargin), 2) == 0, 'params_bus:args', 'Expected name/value pairs.');
for k = 1:2:numel(varargin)
    assert(isfield(p, varargin{k}), 'params_bus:unknown', ...
        'Unknown parameter "%s".', varargin{k});
    p.(varargin{k}) = varargin{k+1};
end

% derived
p.m          = p.m_curb + p.n_pass * p.m_pass;
p.m_eff      = p.m * (1 + p.eps_rot);
p.usable_kwh = p.capacity_kwh * p.soc_window;
end

function F = road_load(v, a, grade, p)
%ROAD_LOAD Tractive force and wheel power from the road load equation.
%   F = ROAD_LOAD(V, A, GRADE, P): V [m/s], A [m/s^2], GRADE as rise/run
%   (scalar or per-sample vector), P from PARAMS_BUS.
%
%       F_roll  = Crr * m * g * cos(theta)
%       F_aero  = 0.5 * rho * Cd * Af * v^2
%       F_grade = m * g * sin(theta)
%       F_acc   = m_eff * a
%       P_wheel = (F_roll + F_aero + F_grade + F_acc) * v
%
%   The components come back separately, which is what makes the energy-balance
%   test possible and what makes the force-vs-speed plot worth looking at.
%
%   Sign: positive = power the road demands, negative = the bus giving it back.

v = v(:);
a = a(:);
if isscalar(grade); grade = repmat(grade, size(v)); else; grade = grade(:); end
assert(numel(grade) == numel(v), 'road_load:grade', 'grade must be scalar or match v.');

theta = atan(grade);

% Rolling resistance opposes motion, so it must vanish at a standstill --
% without sign(v) a parked bus shows a phantom 1.4 kN of drag.
F.roll  = p.Crr * p.m * p.g * cos(theta) .* sign(v);
F.aero  = 0.5 * p.rho * p.Cd * p.A * v.^2 .* sign(v);
F.grade = p.m * p.g * sin(theta) .* ones(size(v));
F.acc   = p.m_eff * a;

F.total   = F.roll + F.aero + F.grade + F.acc;
F.P_wheel = F.total .* v;

F.P_roll  = F.roll  .* v;
F.P_aero  = F.aero  .* v;
F.P_grade = F.grade .* v;
F.P_acc   = F.acc   .* v;
end

function pt = powertrain(P_wheel, v, p)
%POWERTRAIN Wheel power -> battery power.
%   Motoring:  P_batt = P_wheel / eta_dt         + P_aux
%   Braking:   P_batt = -P_recovered * eta_regen + P_aux
%
%   Three separate things stop regen recovering all of the braking energy, and
%   they are tracked separately because they are different problems:
%     1. low-speed cutoff -- below v_regen_cutoff the friction brakes take over
%     2. power cap        -- braking harder than P_regen_max spills to friction
%     3. path efficiency  -- eta_regen on whatever is actually captured
%   "Regen is 85 % efficient" and "regen recovers 85 % of braking energy" are
%   different claims, and this is what lets the model tell them apart.

P_wheel = P_wheel(:);
v = v(:);

motoring = P_wheel > 0;
P_brake  = -P_wheel .* ~motoring;             % >= 0, braking power at the wheel

if p.regen_enabled
    lost_cutoff = P_brake .* (v < p.v_regen_cutoff);
    capturable  = P_brake - lost_cutoff;
    captured    = min(capturable, p.P_regen_max);
    lost_cap    = capturable - captured;
else
    lost_cutoff = P_brake;
    lost_cap    = zeros(size(P_brake));
    captured    = zeros(size(P_brake));
end

pt.P_trac        = (P_wheel .* motoring) / p.eta_dt;
pt.P_regen       = captured * p.eta_regen;
pt.P_lost_eff    = captured * (1 - p.eta_regen);
pt.P_lost_cutoff = lost_cutoff;
pt.P_lost_cap    = lost_cap;
pt.P_brake_wheel = P_brake;
pt.P_aux         = p.P_aux * ones(size(v));
pt.P_batt        = pt.P_trac - pt.P_regen + pt.P_aux;
pt.motoring      = motoring;

% Can the motor actually deliver this trace? Report it rather than silently
% clipping the speed profile the model was told to follow.
pt.P_wheel_peak = max(P_wheel);
pt.n_over_rated = sum(P_wheel > p.P_mot_max);
pt.feasible     = pt.n_over_rated == 0;
end

function res = simulate_cycle(cyc, p)
%SIMULATE_CYCLE Energy and range for one drive cycle.
%   RES = SIMULATE_CYCLE(CYC, P), CYC from LOAD_CYCLE, P from PARAMS_BUS.
%
%   Quasi-static: the speed trace is taken as given and the model asks what
%   power it demands, sample by sample. Integration is trapezoidal.

if nargin < 2; p = params_bus(); end

F  = road_load(cyc.v, cyc.a, p.grade, p);
pt = powertrain(F.P_wheel, cyc.v, p);
t  = cyc.t;
kWh = @(P) trapz(t, P) / 3.6e6;

% where the energy goes
res.E_batt_kwh  = kWh(pt.P_batt);             % net draw from the pack
res.E_trac_kwh  = kWh(pt.P_trac);
res.E_regen_kwh = kWh(pt.P_regen);
res.E_aux_kwh   = kWh(pt.P_aux);

% road-load split, motoring only, at the wheel
res.E_roll_kwh  = kWh(F.P_roll  .* pt.motoring);
res.E_aero_kwh  = kWh(F.P_aero  .* pt.motoring);
res.E_grade_kwh = kWh(F.P_grade .* pt.motoring);
res.E_acc_kwh   = kWh(F.P_acc   .* pt.motoring);
res.E_wheel_kwh = kWh(F.P_wheel .* pt.motoring);

% what braking offered, and what got away
res.E_brake_kwh       = kWh(pt.P_brake_wheel);
res.E_lost_cutoff_kwh = kWh(pt.P_lost_cutoff);
res.E_lost_cap_kwh    = kWh(pt.P_lost_cap);
res.E_lost_eff_kwh    = kWh(pt.P_lost_eff);

% headline numbers
res.distance_km = cyc.distance_km;
res.kwh_per_km  = res.E_batt_kwh / cyc.distance_km;
res.kwh_per_mi  = res.kwh_per_km * 1.609344;

% A steep enough descent makes the cycle net-regenerating; range is then not a
% number, and reporting a negative one would be worse than reporting none.
if res.kwh_per_km > 0
    res.range_km = p.usable_kwh / res.kwh_per_km;
else
    res.range_km = Inf;
end

res.regen_frac_traction = safediv(res.E_regen_kwh, res.E_trac_kwh);
res.recovery_of_braking = safediv(res.E_regen_kwh, res.E_brake_kwh);
res.aux_frac            = safediv(res.E_aux_kwh,   res.E_batt_kwh);

res.P_wheel_peak_kw = pt.P_wheel_peak / 1e3;
res.P_batt_peak_kw  = max(pt.P_batt) / 1e3;
res.P_regen_peak_kw = max(pt.P_regen) / 1e3;
res.feasible        = pt.feasible;
res.n_over_rated    = pt.n_over_rated;

res.cycle  = cyc.name;
res.label  = p.label;
res.params = p;
res.cyc    = cyc;
res.F      = F;
res.pt     = pt;
res.E_batt_trace_kwh = cumtrapz(t, pt.P_batt) / 3.6e6;
end

function r = safediv(a, b)
if b == 0; r = 0; else; r = a / b; end
end

function cyc = make_cycle(t, v_kph, name)
%MAKE_CYCLE Build a validated cycle struct from time and speed vectors.
%   CYC = MAKE_CYCLE(T, V_KPH, NAME)
%
%   Shared by LOAD_CYCLE (real published traces) and by the tests and the
%   steady-cruise case, which need traces built in memory. Keeping the
%   validation and the derived quantities in one place means a synthetic trace
%   is held to exactly the same standard as a downloaded one.

t = t(:);
v_kph = v_kph(:);

assert(numel(t) == numel(v_kph), 'make_cycle:size', 't and v must be the same length.');
assert(~any(isnan([t; v_kph])),  'make_cycle:nan',  '%s: NaNs in trace.', name);
assert(all(diff(t) > 0),         'make_cycle:time', '%s: time is not strictly increasing.', name);
assert(all(v_kph >= 0),          'make_cycle:neg',  '%s: negative speed.', name);
assert(numel(t) >= 10,           'make_cycle:short','%s: trace too short.', name);

dt = diff(t);
if max(abs(dt - dt(1))) > 1e-9
    warning('make_cycle:nonUniform', '%s: non-uniform sample rate (%.3f-%.3f s).', ...
        name, min(dt), max(dt));
end

cyc.name  = name;
cyc.t     = t;
cyc.dt    = dt(1);
cyc.v     = v_kph / 3.6;                       % m/s, SI internally

% Central difference: keeps a the same length as v and, unlike diff(), does
% not shift the acceleration half a sample relative to the speed it came from.
cyc.a = gradient(cyc.v, cyc.t);

cyc.duration_s  = t(end) - t(1);
cyc.distance_km = trapz(t, cyc.v) / 1000;
cyc.v_avg_kph   = cyc.distance_km / (cyc.duration_s/3600);
cyc.v_max_kph   = max(v_kph);
cyc.a_max       = max(cyc.a);
cyc.a_min       = min(cyc.a);

% A stop is a moving -> stationary transition, not every stationary sample.
moving = cyc.v > 0.1;
cyc.n_stops      = sum(diff(moving) == -1);
cyc.stops_per_km = cyc.n_stops / max(cyc.distance_km, eps);
cyc.idle_frac    = mean(~moving);
end

function cyc = load_cycle(name)
%LOAD_CYCLE Build a cycle struct from one of the embedded published traces.
%   CYC = LOAD_CYCLE('manhattan') uses the trace held in CYCLE_MPH10.
%   CYC = LOAD_CYCLE('/path/to/file.csv') still reads a two-column CSV
%   (time_s, speed_kph), so your own traces work unchanged.
%
%   MAKE_CYCLE does the validation and derives distance, stop count and
%   acceleration, so an embedded trace is held to the same standard as a file.

if exist(name, 'file') == 2
    M = readmatrix(name);
    assert(size(M,2) >= 2, 'load_cycle:cols', ...
        '%s: expected columns time_s, speed_kph.', name);
    [~, cycname] = fileparts(name);
    cyc = make_cycle(M(:,1), M(:,2), cycname);
    cyc.file = name;
    return
end

v_kph = double(cycle_mph10(name)) * 0.1 * 1.609344;
t = (0:numel(v_kph)-1)';
cyc = make_cycle(t, v_kph, lower(name));
cyc.file = '<embedded>';
end

function S = sweep(cyc, param, values, varargin)
%SWEEP Re-run one cycle across a range of one parameter.
%   S = SWEEP(CYC, 'n_pass', 0:10:80) sweeps passenger load.
%   S = SWEEP(CYC, 'grade', -0.06:0.01:0.06, 'P_aux', 20e3) sweeps grade with
%   the auxiliary load held at a non-default value.
%
%   Returns column vectors so the result plots directly:
%     S.values .kwh_per_km .kwh_per_mi .range_km .regen_kwh .recovery_of_braking
%     .feasible

n = numel(values);
S.param  = param;
S.values = values(:);
S.cycle  = cyc.name;

fields = {'kwh_per_km','kwh_per_mi','range_km','E_batt_kwh','E_regen_kwh', ...
          'recovery_of_braking','regen_frac_traction','aux_frac','P_batt_peak_kw'};
for f = fields; S.(f{1}) = zeros(n,1); end
S.feasible = true(n,1);

for k = 1:n
    v = values(k);
    if iscell(values); v = values{k}; end
    p = params_bus(varargin{:}, param, v);
    r = simulate_cycle(cyc, p);
    for f = fields; S.(f{1})(k) = r.(f{1}); end
    S.feasible(k) = r.feasible;
end

% Quote everything against the UN-SWEPT baseline, not against the first point
% of the sweep. For a grade sweep running from -6 % to +6 % the first point is
% a downhill run, and quoting the rest against it produces nonsense.
base = simulate_cycle(cyc, params_bus(varargin{:}));
S.base_kwh_per_km = base.kwh_per_km;
S.base_range_km   = base.range_km;
S.kwh_per_km_rel  = S.kwh_per_km / base.kwh_per_km - 1;
S.range_rel       = S.range_km   / base.range_km   - 1;
end


% ===========================================================================
% Drive cycle data
% ===========================================================================

function v = cycle_mph10(name)
%CYCLE_MPH10 Published drive-cycle speed traces, as integer tenths of a mph.
%   The five cycles are embedded rather than shipped as CSVs so the project
%   stays two files. Every trace is uniformly sampled at 1 Hz from t = 0, so
%   only the speed vector is needed; LOAD_CYCLE rebuilds the time vector.
%
%   Stored as tenths of a mph because that is the resolution the sources are
%   published at -- the round trip back to km/h is exact to 4e-07 km/h, far
%   below the precision of anything downstream.
%
%   Provenance, unchanged from the CSVs these replace:
%     manhattan  Manhattan Bus Cycle, NREL DriveCAT (CYC_Manhattan)
%     nybus      New York Bus Cycle, NREL DriveCAT (CYC_NewYorkBus)
%     cbd        Central Business District cycle, NREL DriveCAT (CYC_CBDBUS)
%     udds       EPA Urban Dynamometer Driving Schedule (uddscol.txt)
%     hwfet      EPA Highway Fuel Economy Test (hwycol.txt)

switch lower(name)
    case 'cbd'
        v = [
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 25 48 70 90 111 130 149 167 184 200 200 200 ...
    200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 162 118 73 29 0 0 0 0 0 ...
    0 0 11 34 57 78 99 119 138 156 174 191 200 200 200 200 200 200 200 200 200 200 200 200 ...
    200 200 200 200 200 200 189 145 100 56 11 0 0 0 0 0 0 0 21 43 65 86 107 126 145 163 181 ...
    198 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 171 127 82 ...
    38 0 0 0 0 0 0 0 7 30 52 74 95 115 134 153 170 188 200 200 200 200 200 200 200 200 200 ...
    200 200 200 200 200 200 200 200 200 200 153 109 65 20 0 0 0 0 0 0 0 16 39 61 82 103 122 ...
    141 160 177 195 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 ...
    180 136 91 47 2 0 0 0 0 0 0 2 25 48 70 90 111 130 149 167 184 200 200 200 200 200 200 ...
    200 200 200 200 200 200 200 200 200 200 200 200 200 162 118 73 29 0 0 0 0 0 0 0 11 34 57 ...
    78 99 119 138 156 174 191 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 ...
    200 200 200 189 145 100 56 11 0 0 0 0 0 0 0 21 43 65 86 107 126 145 163 181 198 200 200 ...
    200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 171 127 82 38 0 0 0 0 0 ...
    0 0 7 30 52 74 95 115 134 153 170 188 200 200 200 200 200 200 200 200 200 200 200 200 ...
    200 200 200 200 200 200 200 153 109 65 20 0 0 0 0 0 0 0 16 39 61 82 103 122 141 160 177 ...
    195 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 180 136 91 ...
    47 2 0 0 0 0 0 0 2 25 48 70 90 111 130 149 167 184 200 200 200 200 200 200 200 200 200 ...
    200 200 200 200 200 200 200 200 200 200 162 118 73 29 0 0 0 0 0 0 0 11 34 57 78 99 119 ...
    138 156 174 191 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 ...
    189 145 100 56 11 0 0 0 0 0 0 0 21 43 65 86 107 126 145 163 181 198 200 200 200 200 200 ...
    200 200 200 200 200 200 200 200 200 200 200 200 200 171 127 82 38 0 0 0 0 0 0 0 7 30 52 ...
    74 95 115 134 153 170 188 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 ...
    200 200 200 200 153 109 65 20 0 0 0 0 0 0 0];
    case 'hwfet'
        v = [
    0 0 0 20 49 81 113 145 173 196 218 240 258 271 280 290 300 307 315 322 329 335 341 346 ...
    349 351 357 359 358 353 349 345 346 348 351 357 361 362 365 367 369 370 370 370 370 370 ...
    370 371 373 378 386 393 400 407 414 422 429 435 440 443 445 448 449 450 451 454 457 460 ...
    463 465 468 469 470 471 472 473 472 471 470 469 469 469 470 471 471 472 471 470 469 465 ...
    463 462 463 465 469 471 474 477 480 482 485 488 491 492 491 491 490 490 491 492 493 494 ...
    495 495 495 494 491 489 486 484 481 477 474 473 475 478 479 480 479 479 479 480 480 480 ...
    479 473 460 433 412 395 392 390 390 391 395 401 410 420 431 437 441 443 444 446 447 449 ...
    452 457 459 463 468 469 470 471 476 479 480 480 479 478 473 467 462 459 457 455 454 453 ...
    450 440 431 422 415 415 421 429 435 439 436 433 430 431 434 439 443 446 449 448 444 439 ...
    434 432 432 431 430 430 431 434 439 440 435 426 415 407 400 400 403 410 420 427 431 432 ...
    434 439 443 447 451 454 458 465 469 472 474 473 473 472 472 472 471 470 470 469 468 469 ...
    470 472 475 479 480 480 480 480 480 481 482 482 481 486 489 491 491 491 491 491 490 489 ...
    482 477 475 472 467 462 460 458 456 454 452 450 447 445 442 435 428 420 401 386 375 358 ...
    347 340 333 325 317 306 296 288 284 286 295 314 334 356 375 391 402 411 418 424 428 433 ...
    438 443 447 450 452 454 455 458 460 461 465 468 471 477 483 490 497 503 510 517 524 531 ...
    538 545 552 558 564 569 570 571 573 576 578 580 581 584 587 588 589 590 590 589 588 586 ...
    584 582 581 580 579 576 574 572 571 570 570 569 569 569 570 570 570 570 570 570 570 570 ...
    570 569 568 565 562 560 560 560 561 564 567 569 571 573 574 574 572 570 569 566 563 561 ...
    564 567 571 575 578 580 580 580 580 580 580 579 578 577 577 578 579 580 581 584 589 591 ...
    594 598 599 599 598 596 594 592 591 590 589 587 586 585 584 584 583 582 581 580 579 579 ...
    579 579 579 580 581 581 582 582 582 581 580 580 580 580 580 580 579 579 580 581 581 582 ...
    583 583 583 582 581 580 578 575 571 570 566 561 560 558 555 552 551 550 549 549 549 549 ...
    549 549 550 550 550 550 550 550 551 551 550 549 549 548 547 546 544 543 543 542 541 541 ...
    541 540 540 540 540 540 540 540 540 541 542 545 548 549 550 551 552 552 553 554 555 556 ...
    557 558 559 560 560 560 560 560 560 560 560 560 560 560 560 560 560 559 559 559 558 556 ...
    554 552 551 550 549 546 544 542 541 538 534 533 531 529 526 524 522 521 520 520 520 520 ...
    521 520 520 519 516 514 511 507 503 498 493 487 482 481 480 480 481 484 489 490 491 491 ...
    490 490 489 486 483 480 479 478 477 479 483 490 491 490 489 480 471 462 461 461 462 469 ...
    478 490 497 506 515 522 527 530 536 540 541 544 547 551 554 554 550 545 536 525 502 482 ...
    465 462 460 460 463 468 475 482 488 495 502 507 511 517 522 525 521 516 511 510 510 511 ...
    514 517 520 522 525 528 527 526 523 523 524 525 527 527 524 521 517 511 505 501 498 497 ...
    496 495 495 497 500 502 506 511 516 519 520 521 524 529 533 537 542 545 548 550 555 559 ...
    561 563 564 565 567 569 570 573 577 582 588 591 592 591 588 585 581 577 573 571 568 565 ...
    562 555 546 541 537 532 529 525 520 513 505 495 485 476 468 456 442 425 392 359 326 293 ...
    268 245 215 195 174 151 124 97 70 50 33 20 7 0 0 0];
    case 'manhattan'
        v = [
    0 0 0 0 0 0 0 0 0 0 0 3 6 28 49 71 117 134 149 177 201 206 209 214 235 240 243 239 235 ...
    220 194 171 147 101 64 47 33 22 11 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 17 40 59 73 ...
    77 81 88 95 112 113 115 109 104 110 115 119 125 123 115 113 115 119 112 110 111 118 127 ...
    131 124 97 62 35 11 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16 35 54 81 110 143 167 ...
    180 191 211 228 234 238 224 207 189 162 146 131 99 63 42 22 11 1 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 15 33 49 68 86 103 116 118 118 116 113 115 116 118 91 62 41 29 36 46 55 ...
    45 22 15 10 18 31 47 59 84 96 108 110 107 106 105 103 101 102 93 78 67 72 82 95 108 119 ...
    127 134 141 132 134 130 108 77 44 18 21 28 22 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    13 27 68 95 119 159 193 210 225 242 242 228 213 196 179 141 91 68 48 25 2 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 16 27 30 38 51 70 84 99 113 130 145 155 166 175 171 179 168 164 ...
    161 164 164 161 160 147 131 117 94 68 41 20 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 23 ...
    33 41 59 86 113 155 170 182 201 228 241 253 250 238 225 212 193 174 140 96 68 42 17 1 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 24 48 75 100 119 132 141 138 133 127 116 95 72 41 21 ...
    13 11 15 23 30 45 62 75 84 83 84 86 88 86 70 65 65 68 54 38 41 39 40 43 42 24 8 4 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10 26 42 74 105 124 141 153 164 189 211 223 233 242 233 ...
    221 202 182 155 127 71 22 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 13 24 36 54 67 79 90 100 ...
    107 115 118 125 129 137 141 145 148 150 153 155 157 157 152 156 155 156 154 147 152 153 ...
    153 153 153 153 152 153 154 151 155 160 164 170 168 174 179 168 154 136 112 83 49 16 1 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 6 28 49 71 117 134 149 177 201 206 209 214 235 240 ...
    243 239 235 220 194 171 147 101 64 47 33 22 11 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    17 40 59 73 77 81 88 95 112 113 115 109 104 110 115 119 125 123 115 113 115 119 112 110 ...
    111 118 127 131 124 97 62 35 11 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16 35 54 81 ...
    110 143 167 180 191 211 228 234 238 224 207 189 162 146 131 99 63 42 22 11 1 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 15 33 49 68 86 103 116 118 118 116 113 115 116 118 91 62 41 ...
    29 36 46 55 45 22 15 10 18 31 47 59 84 96 108 110 107 106 105 103 101 102 93 78 67 72 82 ...
    95 108 119 127 134 141 132 134 130 108 77 44 18 21 28 22 8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 13 27 68 95 119 159 193 210 225 242 242 228 213 196 179 141 91 68 48 25 2 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16 27 30 38 51 70 84 99 113 130 145 155 166 175 171 ...
    179 168 164 161 164 164 161 160 147 131 117 94 68 41 20 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 2 23 33 41 59 86 113 155 170 182 201 228 241 253 250 238 225 212 193 174 140 96 ...
    68 42 17 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 24 48 75 100 119 132 141 138 133 127 ...
    116 95 72 41 21 13 11 15 23 30 45 62 75 84 83 84 86 88 86 70 65 65 68 54 38 41 39 40 43 ...
    42 24 8 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10 26 42 74 105 124 141 153 164 189 211 ...
    223 233 242 233 221 202 182 155 127 71 22 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 13 24 36 ...
    54 67 79 90 100 107 115 118 125 129 137 141 145 148 150 153 155 157 157 152 156 155 156 ...
    154 147 152 153 153 153 153 153 152 153 154 151 155 160 164 170 168 174 179 168 154 136 ...
    112 83 49 16 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    case 'nybus'
        v = [
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 96 108 106 102 96 90 92 108 126 114 78 42 30 18 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 96 150 158 154 142 134 142 162 ...
    158 154 144 136 138 174 212 240 276 300 308 306 282 268 282 296 300 278 260 244 198 160 ...
    134 120 78 50 30 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 96 110 106 98 64 ...
    30 28 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 96 ...
    146 148 144 138 136 138 154 152 142 122 88 66 30 28 8 0 0 0 0 0 14 64 100 102 82 54 46 ...
    44 52 76 92 88 74 42 30 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 48 96 150 152 146 136 138 144 158 144 110 72 30 28 6 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 96 158 164 156 144 134 136 162 160 152 144 134 ...
    108 68 30 28 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 48 96 108 106 100 94 88 84 86 118 134 120 84 42 30 26 6 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 96 110 106 98 64 30 28 10 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 96 136 138 134 124 138 154 152 142 ...
    138 146 174 196 194 188 182 200 198 192 184 176 162 126 96 62 30 28 22 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 96 110 108 ...
    102 96 92 94 114 130 124 78 42 30 24 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0];
    case 'udds'
        v = [
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 30 59 86 115 143 169 173 181 207 217 224 225 ...
    221 215 209 204 198 170 149 149 152 155 160 171 191 211 227 229 227 226 213 190 171 158 ...
    158 177 198 216 232 242 246 249 250 246 245 247 248 247 246 246 251 256 257 254 249 250 ...
    254 260 260 257 261 267 275 286 293 298 301 304 307 307 305 304 303 304 308 304 299 295 ...
    298 303 307 309 310 309 304 298 299 302 307 312 318 322 324 322 317 286 253 220 187 154 ...
    121 88 55 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 33 66 99 132 165 198 222 243 258 264 257 251 247 250 252 254 258 272 265 240 227 194 ...
    177 172 181 186 200 222 245 273 305 335 362 373 393 405 421 435 451 460 468 475 475 473 ...
    472 470 470 470 470 470 472 474 479 485 491 495 500 506 510 515 522 532 541 546 549 550 ...
    549 546 546 548 551 555 557 561 563 566 567 567 565 565 565 565 565 565 564 561 558 551 ...
    546 542 540 537 536 539 540 541 541 538 534 530 526 521 524 520 519 517 515 516 518 521 ...
    525 530 535 540 549 554 556 560 560 558 552 545 536 525 515 515 515 511 501 500 501 500 ...
    496 495 495 495 491 486 481 472 461 450 438 426 415 403 385 370 352 338 325 315 306 305 ...
    300 290 275 248 215 201 191 185 170 155 125 108 80 47 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10 ...
    43 76 109 142 173 200 225 237 252 266 281 300 308 316 321 328 336 345 346 349 348 345 ...
    347 355 360 360 360 360 360 360 361 364 365 364 360 351 341 335 314 290 257 230 203 175 ...
    145 120 87 54 21 0 0 0 0 0 0 26 59 92 125 158 191 224 250 256 275 290 300 301 300 297 ...
    293 288 280 250 217 184 151 118 85 52 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 33 66 99 ...
    132 165 198 231 264 278 291 315 330 336 348 351 356 361 360 361 362 360 357 360 360 356 ...
    355 354 352 352 352 352 352 352 350 351 352 355 352 350 350 350 348 346 345 335 320 301 ...
    280 255 225 198 165 132 103 72 40 10 0 0 0 0 0 0 12 35 55 65 85 96 105 119 140 160 177 ...
    190 201 210 220 230 238 245 249 250 250 250 250 250 250 256 258 260 256 252 250 250 250 ...
    244 231 198 165 132 99 66 33 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 33 66 99 130 146 160 170 ...
    170 170 175 177 177 175 170 169 166 170 171 170 166 165 165 166 170 176 185 192 202 210 ...
    211 212 216 220 224 225 225 225 227 237 251 260 265 270 261 228 195 162 129 96 63 30 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 45 78 102 125 140 153 175 196 210 222 ...
    233 245 253 256 260 261 262 262 264 265 265 260 255 236 214 185 164 145 116 87 58 35 20 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 14 33 44 65 92 113 135 146 164 167 165 165 182 192 201 215 ...
    225 225 221 227 233 235 225 216 205 180 150 120 90 62 45 30 21 5 5 32 65 96 125 140 160 ...
    180 196 215 231 245 255 265 271 276 279 283 286 286 283 282 280 275 268 255 235 215 190 ...
    165 149 125 94 62 30 15 15 5 0 30 63 96 129 158 175 184 195 207 220 232 250 265 275 280 ...
    283 289 289 289 288 285 283 283 283 282 276 275 275 275 275 275 275 276 280 285 300 310 ...
    320 330 330 336 340 343 342 340 340 339 336 331 330 325 320 319 316 315 306 300 299 299 ...
    299 299 296 295 295 293 289 282 277 270 255 237 220 205 192 192 201 209 214 220 226 232 ...
    240 250 260 266 266 268 270 272 278 281 288 289 290 291 290 281 275 270 258 250 245 248 ...
    251 255 257 262 269 275 278 284 290 292 291 290 289 285 281 280 280 276 272 266 270 275 ...
    278 280 278 280 280 280 277 274 269 266 265 265 265 263 262 262 259 256 256 259 258 255 ...
    246 235 222 216 216 217 226 234 240 242 244 249 251 252 253 255 252 250 250 250 247 245 ...
    243 243 245 250 250 246 246 241 245 251 256 251 240 220 201 169 136 103 70 37 4 0 0 0 20 ...
    53 86 119 152 175 186 200 211 220 230 245 263 275 281 284 285 285 285 277 275 272 268 ...
    265 260 257 252 240 220 215 215 218 225 230 228 228 230 227 227 227 235 240 246 248 251 ...
    255 256 255 250 241 237 232 229 225 220 216 205 175 142 109 76 43 10 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 12 40 73 106 139 170 185 200 218 230 240 248 256 ...
    265 268 274 279 283 280 275 270 270 263 245 225 215 206 180 150 123 111 106 100 95 91 87 ...
    86 88 90 87 86 80 70 50 42 26 10 0 1 6 16 36 69 100 128 140 145 160 181 200 210 212 213 ...
    214 217 225 230 238 245 250 249 248 250 254 258 260 264 266 269 270 270 270 269 268 268 ...
    265 264 260 255 246 235 215 200 175 160 140 107 74 41 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    21 54 87 120 153 186 211 230 235 230 225 200 167 134 101 68 35 2 0 0 0 0 0 0 0 0 0 0 2 ...
    15 35 65 98 120 129 130 126 128 131 131 140 155 170 186 197 210 215 218 218 215 212 215 ...
    218 220 219 217 215 215 214 201 195 192 196 198 200 195 175 155 130 100 80 60 40 25 7 0 ...
    0 0 0 0 0 0 0 10 10 10 10 10 16 30 40 50 63 80 100 105 95 85 76 88 110 140 170 195 210 ...
    218 222 230 236 241 245 245 240 235 235 235 235 235 235 240 241 245 247 250 254 256 257 ...
    260 262 270 278 283 290 291 290 280 247 214 181 148 115 82 49 16 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 15 48 81 114 132 151 168 183 195 203 213 219 221 224 220 216 ...
    211 205 200 196 185 175 165 155 140 110 80 52 25 0 0 0];
    otherwise
        error('ebus_range:cycle', ...
              'Unknown cycle "%s". Try: manhattan, nybus, cbd, udds, hwfet.', name);
end
v = v(:);
end


% ===========================================================================
% Plotting
% ===========================================================================

function save_fig(f, outfile)
%SAVE_FIG Write a figure to figures/ and close it.
if nargin < 2 || isempty(outfile); return; end
if ~contains(outfile, filesep)
    d = fullfile(pwd, 'figures');
    if ~exist(d, 'dir'); mkdir(d); end
    outfile = fullfile(d, outfile);
end
exportgraphics(f, outfile, 'Resolution', 140);
close(f);
fprintf('  figure: %s\n', outfile);
end

function plot_cycle(res, outfile)
%PLOT_CYCLE Speed, pack power and energy on one shared time axis.
%   The shared axis is the point: a hard stop in the speed trace, the regen
%   spike it produces and the step it takes out of the energy curve all line up.

t = res.cyc.t;
f = figure('Visible', 'off', 'Position', [100 100 1000 780]);
tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, t, res.cyc.v*3.6, 'LineWidth', 1.1, 'Color', [0.10 0.35 0.65]);
ylabel(ax1, 'speed (km/h)'); grid(ax1, 'on');
title(ax1, sprintf('%s  |  %.2f km, %d stops, %.1f km/h avg', ...
    upper(res.cycle), res.distance_km, res.cyc.n_stops, res.cyc.v_avg_kph));

ax2 = nexttile;
P = res.pt.P_batt/1e3;
hold(ax2, 'on');
area(ax2, t, max(P,0), 'FaceColor', [0.80 0.33 0.20], 'EdgeColor', 'none', 'FaceAlpha', 0.75);
area(ax2, t, min(P,0), 'FaceColor', [0.20 0.60 0.35], 'EdgeColor', 'none', 'FaceAlpha', 0.75);
yline(ax2, 0, 'k-'); hold(ax2, 'off');
ylabel(ax2, 'pack power (kW)'); grid(ax2, 'on');
legend(ax2, {'discharge', 'regen'}, 'Location', 'northeast', 'Box', 'off');

% Third panel: the same run with regen switched off. The gap between the two
% curves IS the recovered energy, and it grows at every stop.
p_off = res.params; p_off.regen_enabled = false;
res_off = simulate_cycle(res.cyc, p_off);

ax3 = nexttile;
hold(ax3, 'on');
fill(ax3, [t; flipud(t)], [res.E_batt_trace_kwh; flipud(res_off.E_batt_trace_kwh)], ...
     [0.20 0.60 0.35], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
plot(ax3, t, res_off.E_batt_trace_kwh, '--', 'LineWidth', 1.2, 'Color', [0.45 0.45 0.45]);
plot(ax3, t, res.E_batt_trace_kwh, '-', 'LineWidth', 1.5, 'Color', [0.10 0.35 0.65]);
hold(ax3, 'off');
ylabel(ax3, 'energy drawn (kWh)'); xlabel(ax3, 'time (s)'); grid(ax3, 'on');
legend(ax3, {sprintf('recovered: %.2f kWh (%.0f %% of traction)', ...
             res.E_regen_kwh, 100*res.regen_frac_traction), ...
             sprintf('without regen: %.2f kWh', res_off.E_batt_kwh), ...
             sprintf('with regen: %.2f kWh', res.E_batt_kwh)}, ...
       'Location', 'northwest', 'Box', 'off');
title(ax3, sprintf('%.3f kWh/km  (%.2f kWh/mi)   ->   %.0f km range on %.0f kWh usable', ...
    res.kwh_per_km, res.kwh_per_mi, res.range_km, res.params.usable_kwh));

linkaxes([ax1 ax2 ax3], 'x'); xlim(ax1, [t(1) t(end)]);
save_fig(f, outfile);
end

function plot_forces(res, outfile)
%PLOT_FORCES Where the tractive effort goes.
%   Left: the force components over the cycle. Right: the steady-state
%   force-vs-speed curves. Which side of the crossover in the right-hand panel
%   a route sits on decides everything about its energy signature.

p = res.params; t = res.cyc.t; F = res.F;
f = figure('Visible', 'off', 'Position', [100 100 1150 460]);
tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, t, F.acc/1e3, 'LineWidth', 0.9); hold(ax1, 'on');
plot(ax1, t, F.roll/1e3, 'LineWidth', 1.2);
plot(ax1, t, F.aero/1e3, 'LineWidth', 1.2);
yline(ax1, 0, 'k-'); hold(ax1, 'off'); grid(ax1, 'on');
xlim(ax1, [t(1) t(end)]);
xlabel(ax1, 'time (s)'); ylabel(ax1, 'force (kN)');
title(ax1, sprintf('%s: road load components', upper(res.cycle)));
legend(ax1, {'acceleration', 'rolling', 'aero'}, 'Location', 'best', 'Box', 'off');

ax2 = nexttile;
v = linspace(0, 110/3.6, 300);
F_roll = p.Crr*p.m*p.g * ones(size(v));
F_aero = 0.5*p.rho*p.Cd*p.A*v.^2;
plot(ax2, v*3.6, F_roll/1e3, 'LineWidth', 1.4); hold(ax2, 'on');
plot(ax2, v*3.6, F_aero/1e3, 'LineWidth', 1.4);
plot(ax2, v*3.6, (F_roll+F_aero)/1e3, 'k--', 'LineWidth', 1.2);

% Crossover speed, where aero drag overtakes rolling resistance:
%   v_x = sqrt(2*Crr*m*g / (rho*Cd*A))
v_x = sqrt(2*p.Crr*p.m*p.g/(p.rho*p.Cd*p.A));
xline(ax2, v_x*3.6, ':', sprintf('aero = rolling at %.0f km/h', v_x*3.6), ...
      'LineWidth', 1.2, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');
xline(ax2, res.cyc.v_avg_kph, '-.', sprintf('cycle avg %.0f km/h', res.cyc.v_avg_kph), ...
      'LineWidth', 1.0, 'Color', [0.4 0.4 0.4], 'LabelVerticalAlignment', 'bottom', ...
      'LabelHorizontalAlignment', 'right');
hold(ax2, 'off'); grid(ax2, 'on'); xlim(ax2, [0 110]);
xlabel(ax2, 'speed (km/h)'); ylabel(ax2, 'force (kN)');
title(ax2, 'steady-state road load vs speed');
legend(ax2, {'rolling', 'aero', 'total'}, 'Location', 'northwest', 'Box', 'off');

save_fig(f, outfile);
end

function plot_sweep(S, xlab, outfile, xscale)
%PLOT_SWEEP Consumption and range against one swept parameter.
if nargin < 4; xscale = 1; end

f = figure('Visible', 'off', 'Position', [100 100 900 380]);
tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
x = S.values * xscale;

ax1 = nexttile;
plot(ax1, x, S.kwh_per_km, '-o', 'LineWidth', 1.4, 'MarkerSize', 4);
grid(ax1, 'on'); xlabel(ax1, xlab); ylabel(ax1, 'kWh/km');
title(ax1, sprintf('%s: consumption', upper(S.cycle)));

ax2 = nexttile;
rng = S.range_km;
rng(~isfinite(rng)) = NaN;      % net-regenerating points have no finite range
plot(ax2, x, rng, '-o', 'LineWidth', 1.4, 'MarkerSize', 4);
grid(ax2, 'on'); xlabel(ax2, xlab); ylabel(ax2, 'range (km)');
title(ax2, sprintf('range: %.0f to %.0f km', min(rng), max(rng)));

save_fig(f, outfile);
end

function out = export_bms_profile(res, outfile)
%EXPORT_BMS_PROFILE Write a pack-current trace for the bms-sim project.
%   OUT = EXPORT_BMS_PROFILE(RES) writes data/profiles/<cycle>_<label>.csv.
%   OUT = EXPORT_BMS_PROFILE(RES, FILE) writes FILE.
%
%   bms-sim consumes a (time_s, current_A) drive-cycle profile. This model
%   produces battery POWER, so the conversion needs a pack voltage.
%
%   THE ASSUMPTION, STATED PLAINLY: pack voltage is held at its nominal value,
%   so I = P / V_nom. A real pack sags under load and rises under regen, which
%   means this trace understates current at high discharge and overstates it
%   during regen -- by roughly the same fraction as the voltage error, order
%   5-10 % at the extremes. That is acceptable here precisely because bms-sim
%   models the voltage response itself: feeding it a current derived from a
%   voltage it is about to compute would be circular. The honest split is
%   "this project owns the road load, that project owns the electrochemistry".
%
%   Sign convention matches bms-sim: positive current = discharge.

p = res.params;
t = res.cyc.t;
I = res.pt.P_batt / p.V_pack_nom;             % A, + = discharge

if nargin < 2
    d = fullfile(pwd, 'profiles');
    if ~exist(d, 'dir'); mkdir(d); end
    label = regexprep(lower(res.label), '[^a-z0-9]+', '_');
    outfile = fullfile(d, sprintf('%s_%s.csv', res.cycle, label));
end

fid = fopen(outfile, 'w');
fprintf(fid, '# ebus-range pack current profile\n');
fprintf(fid, '# cycle=%s  params=%s\n', res.cycle, res.label);
fprintf(fid, '# V_pack_nom=%.1f V (held constant -- see EXPORT_BMS_PROFILE)\n', p.V_pack_nom);
fprintf(fid, '# net energy=%.3f kWh over %.3f km (%.3f kWh/km)\n', ...
        res.E_batt_kwh, res.distance_km, res.kwh_per_km);
fprintf(fid, '# peak discharge=%.1f A  peak charge=%.1f A  (+ = discharge)\n', max(I), -min(I));
fprintf(fid, 'time_s,current_A\n');
fprintf(fid, '%g,%.4f\n', [t, I]');
fclose(fid);

fprintf('wrote %s  (%d samples, %.1f to %.1f A)\n', outfile, numel(t), min(I), max(I));
out = outfile;
end


% ===========================================================================
% Entry points
% ===========================================================================

function results = run_all()
%RUN_ALL Every cycle, one parameter set, one table.
%   Same bus, same parameters, five published cycles. The spread is the point.


cycles = {'manhattan', 'nybus', 'cbd', 'udds', 'hwfet'};
p = params_bus();
results = cell(1, numel(cycles));

fprintf('\nBus: %.0f kg (%d pax), Cd %.2f, A %.1f m2, Crr %.3f, aux %.1f kW, %.0f kWh usable\n\n', ...
    p.m, p.n_pass, p.Cd, p.A, p.Crr, p.P_aux/1e3, p.usable_kwh);
fprintf('%-11s %7s %7s %7s %9s %9s %8s %8s %7s\n', ...
    'cycle', 'km', 'km/h', 'stop/km', 'kWh/km', 'kWh/mi', 'regen %', 'aux %', 'range');
fprintf('%s\n', repmat('-', 1, 82));

for k = 1:numel(cycles)
    cyc = load_cycle(cycles{k});
    r = simulate_cycle(cyc, p);
    results{k} = r;
    if r.feasible; flag = ''; else; flag = '  [!]'; end
    fprintf('%-11s %7.2f %7.1f %7.1f %9.3f %9.2f %8.0f %8.0f %6.0f km%s\n', ...
        r.cycle, r.distance_km, cyc.v_avg_kph, cyc.stops_per_km, ...
        r.kwh_per_km, r.kwh_per_mi, 100*r.regen_frac_traction, 100*r.aux_frac, ...
        r.range_km, flag);
end
fprintf('%s\n', repmat('-', 1, 82));
fprintf(['  regen %% = energy returned to the pack as a fraction of traction energy drawn\n' ...
         '  [!]     = trace demands more than the %.0f kW rated traction power; see README\n\n'], ...
         p.P_mot_max/1e3);

report(results{1});          % manhattan, in detail
report(results{5});          % hwfet, for contrast

plot_comparison(results, 'cycle_comparison.png');
plot_cycle(results{1},  'cycle_manhattan.png');
plot_cycle(results{5},  'cycle_hwfet.png');
plot_forces(results{1}, 'forces_manhattan.png');
plot_forces(results{5}, 'forces_hwfet.png');

if nargout == 0; clear results; end
end

% ---------------------------------------------------------------------------

function report(r)
fprintf('=== %s ===\n', upper(r.cycle));
fprintf('  consumption  %.3f kWh/km (%.2f kWh/mi),  range %.0f km\n', ...
    r.kwh_per_km, r.kwh_per_mi, r.range_km);
fprintf('  energy       net %.2f kWh = traction %.2f - regen %.2f + aux %.2f\n', ...
    r.E_batt_kwh, r.E_trac_kwh, r.E_regen_kwh, r.E_aux_kwh);
fprintf('  road load    roll %.2f | aero %.2f | accel %.2f kWh at the wheel\n', ...
    r.E_roll_kwh, r.E_aero_kwh, r.E_acc_kwh);
fprintf('  braking      %.2f kWh offered -> %.2f kWh recovered (%.0f %%)\n', ...
    r.E_brake_kwh, r.E_regen_kwh, 100*r.recovery_of_braking);
fprintf('               lost: %.2f cutoff | %.2f power cap | %.2f path efficiency\n', ...
    r.E_lost_cutoff_kwh, r.E_lost_cap_kwh, r.E_lost_eff_kwh);
fprintf('  peaks        wheel %.0f kW | pack %.0f kW | regen %.0f kW\n\n', ...
    r.P_wheel_peak_kw, r.P_batt_peak_kw, r.P_regen_peak_kw);
end

function plot_comparison(results, outfile)
% Built at the PACK so the stack sums exactly to the net consumption beside it:
%   net = (E_roll + E_aero + E_grade + E_acc)/eta_dt - E_regen + E_aux
n = numel(results);
names = cell(1,n); E = zeros(n,4); kwh_km = zeros(1,n);
for k = 1:n
    r = results{k}; d = r.distance_km; e = r.params.eta_dt;
    names{k} = upper(r.cycle);
    E(k,:) = [r.E_roll_kwh/e, r.E_aero_kwh/e, ...
              (r.E_acc_kwh + r.E_grade_kwh)/e - r.E_regen_kwh, r.E_aux_kwh] / d;
    kwh_km(k) = r.kwh_per_km;
    assert(abs(sum(E(k,:)) - kwh_km(k)) < 1e-9*kwh_km(k), ...
        'run_all:balance', '%s: stack does not sum to net consumption.', r.cycle);
end

f = figure('Visible', 'off', 'Position', [100 100 980 400]);
tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
b = bar(ax1, categorical(names, names), E, 'stacked');
b(1).FaceColor = [0.35 0.45 0.60]; b(2).FaceColor = [0.55 0.70 0.85];
b(3).FaceColor = [0.80 0.45 0.30]; b(4).FaceColor = [0.60 0.60 0.60];
ylabel(ax1, 'kWh/km at the pack'); grid(ax1, 'on');
legend(ax1, {'rolling', 'aero', 'accel net of regen', 'auxiliaries'}, ...
       'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
title(ax1, 'where the energy goes (incl. drivetrain loss)');

ax2 = nexttile;
bar(ax2, categorical(names, names), kwh_km, 'FaceColor', [0.20 0.40 0.55]);
hold(ax2, 'on');
yline(ax2, 1.5/1.609344, '--', '1.5 kWh/mi', 'LineWidth', 1.1, 'Color', [0.3 0.5 0.3], ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
yline(ax2, 3.0/1.609344, '--', '3.0 kWh/mi', 'LineWidth', 1.1, 'Color', [0.3 0.5 0.3], ...
      'LabelHorizontalAlignment', 'left');
hold(ax2, 'off');
ylabel(ax2, 'kWh/km at the pack'); grid(ax2, 'on');
title(ax2, 'net consumption vs published transit band');

save_fig(f, outfile);
end

function run_sweeps()
%RUN_SWEEPS Sensitivity of consumption and range to the things that vary in
%   service: how full the bus is, what the road does, how good the regen path
%   is, and how fast you cruise. Everything is quoted against the baseline.


cyc = load_cycle('cbd');
fprintf('\nSensitivity on %s (%.2f km, %.1f km/h avg, %d stops)\n', ...
    upper(cyc.name), cyc.distance_km, cyc.v_avg_kph, cyc.n_stops);

% --- passenger load ------------------------------------------------------
S = sweep(cyc, 'n_pass', 0:10:60);
banner('passenger load');
for k = 1:numel(S.values)
    p = params_bus('n_pass', S.values(k));
    row(sprintf('%d pax (%.1f t)', S.values(k), p.m/1000), S, k);
end
plot_sweep(S, 'passengers', 'sweep_passengers.png');

% --- grade ---------------------------------------------------------------
S = sweep(cyc, 'grade', (-6:2:6)/100);
banner('road grade');
for k = 1:numel(S.values)
    row(sprintf('%+.0f %% grade', 100*S.values(k)), S, k);
end
plot_sweep(S, 'grade (%)', 'sweep_grade.png', 100);
grade_asymmetry(S);

% --- regen path efficiency ----------------------------------------------
S = sweep(cyc, 'eta_regen', 0:0.15:0.9);
banner('regen path efficiency');
for k = 1:numel(S.values)
    row(sprintf('eta_regen = %.2f', S.values(k)), S, k);
end
plot_sweep(S, 'regen path efficiency', 'sweep_eta_regen.png');

% --- regen power cap -----------------------------------------------------
S = sweep(cyc, 'P_regen_max', (0:50:300)*1e3);
banner('regen power cap');
for k = 1:numel(S.values)
    row(sprintf('%d kW cap', round(S.values(k)/1e3)), S, k);
end
plot_sweep(S, 'regen power cap (kW)', 'sweep_regen_cap.png', 1e-3);
fprintf(['  Flattens once the cap clears the cycle''s hardest stop -- past that\n' ...
         '  point a bigger inverter buys nothing on this duty cycle.\n']);

% --- cruise speed --------------------------------------------------------
% Not a drive cycle: a constant-speed operating point, which is what an
% intercity coach approximates for hours at a time. Aero goes as v^3, so aero
% energy per km goes as v^2, and that is the whole highway story.
fprintf('\ncruise speed (steady, flat, no stops)\n');
fprintf('%-22s %9s %9s %9s\n%s\n', 'setting', 'kWh/km', 'range km', 'aero %', repmat('-', 1, 52));
p = params_bus();
speeds = 50:10:110;
kwh = zeros(size(speeds)); rng_km = zeros(size(speeds));
for k = 1:numel(speeds)
    c = make_cycle((0:3600)', speeds(k)*ones(3601,1), sprintf('cruise%d', speeds(k)));
    r = simulate_cycle(c, p);
    kwh(k) = r.kwh_per_km; rng_km(k) = r.range_km;
    fprintf('%-22s %9.3f %9.0f %8.0f%%\n', sprintf('%d km/h', speeds(k)), ...
        kwh(k), rng_km(k), 100*r.E_aero_kwh/(r.E_aero_kwh + r.E_roll_kwh));
end
fprintf('%s\n', repmat('-', 1, 52));
fprintf('  Range is %.0f km at 80 km/h and %.0f km at 100. Cruise speed is a range knob.\n\n', ...
    rng_km(speeds==80), rng_km(speeds==100));

f = figure('Visible', 'off', 'Position', [100 100 900 380]);
tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile; plot(ax1, speeds, kwh, '-o', 'LineWidth', 1.4); grid(ax1, 'on');
xlabel(ax1, 'cruise speed (km/h)'); ylabel(ax1, 'kWh/km');
title(ax1, 'steady-cruise consumption');
ax2 = nexttile; plot(ax2, speeds, rng_km, '-o', 'LineWidth', 1.4); grid(ax2, 'on');
xlabel(ax2, 'cruise speed (km/h)'); ylabel(ax2, 'range (km)');
title(ax2, sprintf('range: %.0f km at 50, %.0f km at 110', rng_km(1), rng_km(end)));
save_fig(f, 'sweep_cruise_speed.png');
end

% ---------------------------------------------------------------------------

function banner(what)
fprintf('\n%s\n%-22s %9s %9s %9s %9s\n%s\n', what, ...
    'setting', 'kWh/km', 'vs base', 'range km', 'regen %', repmat('-', 1, 62));
end

function row(label, S, k)
if isfinite(S.range_km(k)); r = sprintf('%9.0f', S.range_km(k));
else;                       r = sprintf('%9s', 'net regen'); end
fprintf('%-22s %9.3f %+8.0f%% %s %9.0f\n', label, S.kwh_per_km(k), ...
    100*S.kwh_per_km_rel(k), r, 100*S.recovery_of_braking(k));
end

function grade_asymmetry(S)
% What a climb costs is not what the matching descent gives back: the energy
% goes up through the drivetrain and comes back through the regen path, and
% both are lossy. Measure it rather than asserting it.
i0 = find(S.values == 0);
fprintf('\n  A climb costs more than the matching descent returns:\n');
for k = 1:numel(S.values)
    g = S.values(k);
    if g <= 0; continue; end
    j = find(abs(S.values + g) < 1e-12);
    if isempty(j); continue; end
    cost = S.kwh_per_km(k)  - S.kwh_per_km(i0);
    gain = S.kwh_per_km(i0) - S.kwh_per_km(j);
    fprintf('    +/-%.0f %%: climb costs %.3f kWh/km, descent returns %.3f -> %.0f %% back\n', ...
        100*g, cost, gain, 100*gain/cost);
end
end

function export_profiles()
%EXPORT_PROFILES Write pack-current traces for the bms-sim project.
%   Produces data/profiles/*.csv as (time_s, current_A), the form bms-sim's
%   CSV drive-cycle profile expects. The voltage assumption behind the
%   power -> current conversion is in export_bms_profile.m and is repeated in
%   the header of every CSV, so nothing downstream has to guess.


fprintf('\n');
for c = {'manhattan', 'cbd'}
    r = simulate_cycle(load_cycle(c{1}), params_bus());
    export_bms_profile(r);
end
fprintf('\nPositive current is discharge, matching bms-sim''s convention.\n\n');
end


% ===========================================================================
% Tests
% ===========================================================================

function run_tests()
%RUN_TESTS Validation suite for the road-load and energy model.
%   Run from anywhere:  run_tests()
%
%   Three kinds of check:
%     ANALYTIC     -- cases with a closed-form answer the model must reproduce.
%     STRUCTURAL   -- conservation identities that hold for any input, and that
%                     catch a flipped sign instantly.
%     PLAUSIBILITY -- against published real-world consumption figures.


global PASSED FAILED
PASSED = 0; FAILED = 0;

fprintf('\n--- analytic ---------------------------------------------------\n');
test_constant_speed();
test_constant_acceleration();
test_aero_cubed();

fprintf('\n--- structural -------------------------------------------------\n');
test_energy_balance();
test_braking_accounting();
test_battery_accounting();
test_regen_off();
test_grade_round_trip();
test_parked_bus();

fprintf('\n--- plausibility -----------------------------------------------\n');
test_real_world_band();
test_cycle_metadata();

fprintf('\n================================================================\n');
fprintf('  %d passed, %d failed\n\n', PASSED, FAILED);
if FAILED > 0; error('run_tests:failed', '%d test(s) failed.', FAILED); end
end

% =====================================================================
% analytic
% =====================================================================

function test_constant_speed()
% At constant speed on the flat there is no acceleration or grade term, so
% wheel power is exactly (F_roll + F_aero)*v -- no integration error to hide in.
p = params_bus('P_aux', 0);
cyc = make_cycle((0:600)', 60*ones(601,1), 'cruise60');
F = road_load(cyc.v, cyc.a, 0, p);
v = 60/3.6;
expect = (p.Crr*p.m*p.g + 0.5*p.rho*p.Cd*p.A*v^2) * v;
check('constant speed: P_wheel matches closed form', ...
      relerr(mean(F.P_wheel(10:end-10)), expect) < 1e-12);
end

function test_constant_acceleration()
% With drag switched off, accelerating from rest must cost exactly the kinetic
% energy of the EFFECTIVE mass -- which also proves the rotating-inertia
% allowance is really in there.
p = params_bus('Crr', 0, 'Cd', 0, 'P_aux', 0, 'eta_dt', 1);
t = (0:0.1:20)'; a = 0.5;
cyc = make_cycle(t, a*t*3.6, 'ramp');
res = simulate_cycle(cyc, p);
vf = a*t(end);
check('constant accel: energy = 1/2 m_eff v^2', ...
      relerr(res.E_batt_kwh, 0.5*p.m_eff*vf^2/3.6e6) < 1e-6);
end

function test_aero_cubed()
% Aero power goes as v^3, so aero energy over a constant-speed run is
% 1/2 rho Cd A v^3 t.
p = params_bus('Crr', 0, 'P_aux', 0, 'eta_dt', 1);
v = 100/3.6;
cyc = make_cycle((0:900)', 100*ones(901,1), 'cruise100');
res = simulate_cycle(cyc, p);
check('aero: energy = 1/2 rho Cd A v^3 t', ...
      relerr(res.E_aero_kwh, 0.5*p.rho*p.Cd*p.A*v^3*900/3.6e6) < 1e-9);
end

% =====================================================================
% structural
% =====================================================================

function test_energy_balance()
% The four road-load components must sum to the total, pointwise and in the
% integral. Cheap, and the first thing that breaks if a sign flips.
cyc = load_cycle('manhattan');
p = params_bus('grade', 0.02);
F = road_load(cyc.v, cyc.a, p.grade, p);
resid = F.P_wheel - (F.P_roll + F.P_aero + F.P_grade + F.P_acc);
check('components sum to total (pointwise)', max(abs(resid)) < 1e-6*max(abs(F.P_wheel)));

r = simulate_cycle(cyc, p);
check('component integrals sum to wheel energy', ...
      abs(r.E_wheel_kwh - (r.E_roll_kwh + r.E_aero_kwh + r.E_grade_kwh + r.E_acc_kwh)) ...
      < 1e-9*abs(r.E_wheel_kwh));
end

function test_braking_accounting()
% Every joule braking offers is recovered, lost to the low-speed cutoff, lost
% to the power cap, or lost in the path. Nothing may go missing.
r = simulate_cycle(load_cycle('nybus'), params_bus());
total = r.E_regen_kwh + r.E_lost_cutoff_kwh + r.E_lost_cap_kwh + r.E_lost_eff_kwh;
check('braking: recovered + three loss terms = energy offered', ...
      relerr(total, r.E_brake_kwh) < 1e-9);
end

function test_battery_accounting()
r = simulate_cycle(load_cycle('cbd'), params_bus());
check('battery: net = traction - regen + aux', ...
      relerr(r.E_trac_kwh - r.E_regen_kwh + r.E_aux_kwh, r.E_batt_kwh) < 1e-9);
end

function test_regen_off()
% Two ways of switching regen off must agree, and regen must actually help.
cyc = load_cycle('manhattan');
r1 = simulate_cycle(cyc, params_bus('regen_enabled', false));
r2 = simulate_cycle(cyc, params_bus('eta_regen', 0));
check('regen off: disabled == zero efficiency', relerr(r1.E_batt_kwh, r2.E_batt_kwh) < 1e-12);
check('regen on costs less than regen off', ...
      simulate_cycle(cyc, params_bus()).E_batt_kwh < r1.E_batt_kwh);
end

function test_grade_round_trip()
% Climb 100 m at constant speed, then descend the same 100 m. The grade term
% must integrate to zero, but the ROUND TRIP must still cost energy, because
% what goes up through the drivetrain comes back through the regen path.
p = params_bus('P_aux', 0);
v = 40/3.6;
n = round(100/(v*0.05));
t = (0:2*n)';
% Symmetric about the crest, so climb and descent carry equal trapezoidal
% weight -- otherwise the test fails on its own construction.
grade = [0.05*ones(n,1); 0; -0.05*ones(n,1)];
cyc = make_cycle(t, 40*ones(size(t)), 'hill');

F = road_load(cyc.v, cyc.a, grade, p);
check('grade round trip: grade energy integrates to zero', ...
      abs(trapz(t, F.P_grade)) < 1e-9*trapz(t, abs(F.P_grade)));

p.grade = grade;  r_hill = simulate_cycle(cyc, p);
p.grade = 0;      r_flat = simulate_cycle(cyc, p);
check('grade round trip: costs energy once regen is lossy', ...
      r_hill.E_batt_kwh > r_flat.E_batt_kwh * 1.01);
end

function test_parked_bus()
% A parked bus draws exactly the auxiliary load -- no phantom rolling drag.
p = params_bus();
cyc = make_cycle((0:600)', zeros(601,1), 'parked');
F = road_load(cyc.v, cyc.a, 0, p);
check('parked: no phantom road load', max(abs(F.P_wheel)) < 1e-12);
pt = powertrain(F.P_wheel, cyc.v, p);
check('parked: draws exactly P_aux', max(abs(pt.P_batt - p.P_aux)) < 1e-12);
end

% =====================================================================
% plausibility
% =====================================================================

function test_real_world_band()
% NREL fleet evaluations put 40 ft battery-electric transit buses at roughly
% 1.5-3.0 kWh/mi depending on duty cycle and load. A transit cycle landing
% outside that band means the parameters are wrong, not the bus.
for c = {'manhattan', 'cbd', 'nybus'}
    r = simulate_cycle(load_cycle(c{1}), params_bus());
    check(sprintf('%-10s %.2f kWh/mi inside the 1.2-3.5 kWh/mi band', c{1}, r.kwh_per_mi), ...
          r.kwh_per_mi > 1.2 && r.kwh_per_mi < 3.5);
end
end

function test_cycle_metadata()
% Guard against a corrupted or regenerated CSV: each trace must still match
% the published cycle definition it claims to be.
expect = {'manhattan', 1089, 3.32; 'cbd', 574, 3.22; 'nybus', 600, 0.99; ...
          'hwfet', 765, 16.51; 'udds', 1369, 11.99};
for k = 1:size(expect,1)
    cyc = load_cycle(expect{k,1});
    check(sprintf('%-10s matches published duration and distance', expect{k,1}), ...
          abs(cyc.duration_s - expect{k,2}) <= 2 && ...
          abs(cyc.distance_km - expect{k,3}) <= 0.03*expect{k,3});
end
end

% =====================================================================

function check(name, cond)
global PASSED FAILED
if cond
    PASSED = PASSED + 1; fprintf('  PASS  %s\n', name);
else
    FAILED = FAILED + 1; fprintf('  FAIL  %s\n', name);
end
end

function e = relerr(a, b)
e = abs(a - b) / max(abs(b), eps);
end
