# ebus-range

**How far does a battery-electric bus actually go on a given route?**

A drive-cycle energy and range model for a 40 ft battery-electric transit bus,
written in plain MATLAB. Feed it a published transit duty cycle and it returns
consumption in kWh/km, how much energy regenerative braking puts back, and an
estimated range.

Base MATLAB only — no Simulink, no toolboxes. Two files: this README and
`ebus_range.m`, about 1,200 lines including the five drive cycles and the tests.

Figures are not stored in the repo; `ebus_range all` and `ebus_range sweeps`
write them to `./figures`. The ones described below are what those commands
produce.

*Manhattan Bus Cycle (`figures/cycle_manhattan.png`). Top: the published speed trace. Middle: pack power, with
regen below the axis. Bottom: energy drawn, compared against the same run with
regen disabled — the shaded gap is the 3.38 kWh braking put back, 55 % of the
traction energy drawn.*

## Results

Same bus, same parameters, five published cycles:

| cycle | km | km/h avg | stops/km | **kWh/km** | **kWh/mi** | regen % | aux % | range |
|---|---|---|---|---|---|---|---|---|
| Manhattan | 3.32 | 11.0 | 6.0 | 1.274 | 2.05 | 55 | 36 | 311 km |
| NY Bus | 0.99 | 5.9 | 11.1 | 1.827 | 2.94 | 57 | 46 | 217 km ⚠ |
| CBD (SAE J1376) | 3.23 | 20.2 | 4.3 | 0.975 | 1.57 | 49 | 25 | 406 km |
| UDDS | 11.99 | 31.5 | 1.4 | 0.975 | 1.57 | 41 | 16 | 406 km ⚠ |
| HWFET | 16.51 | 77.7 | 0.1 | 1.019 | 1.64 | 12 | 6 | 389 km ⚠ |

17.7 t (40 passengers), Cd 0.60, A 8.0 m², Crr 0.008, 5 kW auxiliaries,
396 kWh usable. `regen %` is energy returned to the pack as a fraction of
traction energy drawn. ⚠ marks cycles that briefly demand more than the 250 kW
rated traction power — see [Limitations](#limitations).

Three things worth pointing at:

- **Regenerative braking matters more in transit duty than people expect.** On
  the Manhattan cycle it returns **55 %** of the traction energy drawn. On the
  highway trace it returns 12 %, because there is nothing to recover from a bus
  that never stops.
- **Auxiliaries are a route-dependent load, not a fixed one.** The same 5 kW
  draw is 6 % of consumption on the highway and 46 % on NY Bus, because
  auxiliaries are drawn against *time* and a transit bus spends most of its time
  barely moving.
- **Transit and highway duty are different problems.** Manhattan spends its
  energy accelerating 17.7 t and gets most of it back. HWFET spends it fighting
  aerodynamic drag and gets almost none of it back.

## How it works

Instantaneous tractive force from the road load equation, integrated over the
cycle:

```
F_roll  = Crr · m · g · cos(θ)          rolling resistance (vanishes at rest)
F_aero  = ½ · ρ · Cd · A · v²           aerodynamic drag
F_grade = m · g · sin(θ)                gravity along the slope
F_acc   = m_eff · a                     m_eff = m·(1+ε), ε = 0.05 rotating inertia

P_wheel = (F_roll + F_aero + F_grade + F_acc) · v
```

Wheel power maps to pack power differently depending on which way it flows:

```
motoring:   P_batt = P_wheel / η_dt          + P_aux
braking:    P_batt = −P_recovered · η_regen  + P_aux
```

Energy is integrated with `cumtrapz`, and acceleration comes from a central
difference (`gradient`), which keeps `a` aligned with the `v` it came from
instead of shifted half a sample.

### Regenerative braking

Three separate things stop regen recovering all of the braking energy, and
`powertrain` tracks them separately because they are
different engineering problems:

1. **Low-speed cutoff** — below ~5 km/h the friction brakes take over. Costs the tail of every stop.
2. **Power cap** — braking harder than the motor's 250 kW rating spills the excess to the friction brakes. Costs the peak of every hard stop.
3. **Path efficiency** — η_regen on whatever is actually captured.

This is what lets the model tell apart two very different claims: *"regen is
85 % efficient"* and *"regen recovers 85 % of braking energy."* On Manhattan it
recovers 83 % of what braking offered; on HWFET, 78 %, because a highway stop
from 90 km/h runs into the power cap.

*`figures/forces_hwfet.png`, right panel: aerodynamic drag overtakes rolling resistance at 78 km/h, where
`v_x = sqrt(2·Crr·m·g / (ρ·Cd·A))`. Which side of that crossover a route sits on
decides its energy signature — HWFET happens to average 77.7 km/h, right on top
of it, while Manhattan averages 11 km/h, so far left that aero never shows up in
its energy budget at all.*

## Validation

Checking the model against reality is the difference between a simulation and a
toy. `ebus_range tests` runs 21 checks in three groups.

**Analytic** — cases with a closed-form answer:
- Constant speed on the flat reproduces `(F_roll + F_aero)·v` to 1e-12.
- Accelerating from rest with drag zeroed costs exactly `½·m_eff·v²`, which also
  proves the rotating-inertia allowance is really in there.
- Aero energy over a constant-speed run equals `½ρCdAv³t`.

**Structural** — identities that must hold for any input:
- The four force components sum to the total, pointwise and in the integral.
- Every joule braking offers is accounted for: recovered, lost to the cutoff,
  lost to the cap, or lost in the path.
- Disabling regen and setting `η_regen = 0` give identical answers.
- A parked bus draws exactly `P_aux` and no phantom rolling drag.
- **Grade round trip**: climb 100 m, descend 100 m. The grade energy integrates
  to zero exactly, but the round trip still costs energy, because what goes up
  through the drivetrain comes back through the lossy regen path.

**Plausibility** — against published data:

> NREL fleet evaluations put 40 ft battery-electric transit buses at roughly
> **1.5–3.0 kWh/mi** depending on duty cycle, load and climate.

The three transit cycles land at **1.57, 2.05 and 2.94 kWh/mi** — inside the
band across its full width, ordered by stop density, which is what you would
expect. If the model said 0.3 or 12, the parameters would be wrong.

The cycle traces themselves are re-checked against their published duration and
distance on every test run, so a corrupted CSV fails loudly rather than quietly
changing every number here.

## Sensitivities

`ebus_range sweeps`, all on the CBD cycle:

| sweep | result |
|---|---|
| passenger load, 0 → 60 pax | +18 % consumption, 454 → 384 km |
| grade, 0 → +6 % | +320 % consumption, 406 → 97 km |
| regen path efficiency, 0 → 0.9 | −44 % consumption, 236 → 424 km |
| regen power cap, 0 → 300 kW | −42 % consumption, flat above 250 kW |
| steady cruise, 50 → 110 km/h | 0.70 → 1.32 kWh/km, 562 → 300 km |

Two results the sweeps make concrete:

- The regen cap curve **flattens once the cap clears the cycle's hardest stop**.
  Past that point a bigger inverter buys nothing on this duty cycle.
- **A climb never fully comes back.** At ±6 % grade, 84 % of what the climb cost
  is returned by the matching descent; at ±2 %, 91 %. The energy goes up through
  the drivetrain and returns through the regen path, and both are lossy.

The cruise sweep is a constant-speed operating point rather than a drive cycle —
which is what an intercity coach approximates for hours at a time. Aero power
goes as v³, so aero energy per km goes as v², and cruise speed becomes a range
knob: 421 km at 80 km/h, 336 km at 100.

## Running it

Nothing to install and nothing to add to the path — `ebus_range.m` is
self-contained. `cd` to the repo and run:

```bash
matlab -batch "ebus_range tests"
```

```bash
matlab -batch "ebus_range all; ebus_range sweeps"
```

Or interactively, `ebus_range` with no argument prints the same list:

| command | what it does |
|---|---|
| `ebus_range all` | every cycle, one parameter set, one table, plus figures |
| `ebus_range sweeps` | passengers, grade, auxiliaries, regen, cruise speed |
| `ebus_range tests` | 21 analytic, structural and plausibility checks |
| `ebus_range profiles` | pack-current CSVs for a battery model |
| `ebus_range help` | the above |

`R = ebus_range('all')` also returns the per-cycle result structs. Figures land
in `./figures` and exported profiles in `./profiles`, both relative to wherever
you run it from.

Inside `ebus_range.m`, in reading order:

```
params_bus           vehicle parameters, every value sourced
road_load            the four force components + wheel power
powertrain           wheel power -> pack power, regen limits
simulate_cycle       integrate; energy, consumption, range
make_cycle           validate a trace, derive distance/stops/acceleration
load_cycle           an embedded cycle, or your own CSV by path
sweep                one-parameter sensitivity driver
cycle_mph10          the five published drive cycles, embedded
save_fig  plot_cycle  plot_forces  plot_sweep  export_bms_profile
run_all  run_sweeps  export_profiles  run_tests
```

## Parameters

| parameter | value | basis |
|---|---|---|
| curb mass | 15,000 kg | typical 40 ft BEB (~33,000 lb; GVWR ~45,000 lb) |
| passengers | 40 × 68 kg | 38 seated + standees; 60 pax ≈ 19.1 t is crush load |
| rotating inertia ε | 0.05 | wheels, gearbox, rotor |
| Cd × A | 0.60 × 8.0 m² | blunt bus body; 2.6 m wide × 3.2 m tall |
| Crr | 0.008 | radial truck tyres, dry pavement |
| η_dt | 0.90 | inverter × motor × gearbox |
| η_regen | 0.85 | wheel → pack on the regen path |
| motor rating | 250 kW | traction and regen |
| auxiliaries | 5 kW | HVAC, air compressor, steering, lighting |
| pack | 440 kWh, 90 % usable | fleets do not run 100–0 % |

These are representative figures from transit-electrification literature, not
any OEM's data. They live in one function, `params_bus`, each
with its basis written next to it — a range number is only as arguable as the
parameters behind it.

## Data

Every trace is a **published** cycle; none were invented or hand-drawn for this
project. All five are embedded in `cycle_mph10`, stored as integer tenths of a
mph — the resolution the sources are published at — and rebuilt to km/h on load.
The round trip is exact to 4e-7 km/h, and every derived quantity (distance, stop
count, idle fraction) is identical to the CSVs they replace.

| cycle | what it is | duration | distance | max speed | source |
|---|---|---|---|---|---|
| `manhattan` | Manhattan Bus Cycle — NYC transit duty, severe stop-and-go | 1089 s | 3.32 km | 40.7 km/h | ADVISOR `CYC_Manhattan` |
| `nybus` | New York Bus Cycle — the most extreme transit cycle in common use | 600 s | 0.99 km | 49.6 km/h | ADVISOR `CYC_NewYorkBus` |
| `cbd` | Central Business District (SAE J1376) — the standard transit cycle | 574 s | 3.23 km | 32.2 km/h | ADVISOR `CYC_CBDBUS` |
| `udds` | Urban Dynamometer Driving Schedule | 1369 s | 11.99 km | 91.3 km/h | EPA `uddscol.txt` |
| `hwfet` | Highway Fuel Economy Test — the highway/coach speed profile | 765 s | 16.51 km | 96.4 km/h | EPA `hwycol.txt` |

Each trace was accepted only if it matched its published duration, maximum
speed and distance to within 3 %. That check is the provenance argument: a trace
reproducing all three headline figures of a published cycle is that cycle,
whatever route it took to get here. It is repeated as a test (`ebus_range
tests`, "matches published duration and distance"), so a corrupted trace fails
loudly rather than quietly changing every result.

The EPA files are authoritative, from
`epa.gov/vehicle-and-fuel-emissions-testing/dynamometer-drive-schedules`. The
three bus cycles come from a redistributed copy of NREL's ADVISOR library rather
than NREL directly (DriveCAT was unreachable at the time of writing), which is
what the 3 % check exists to cover.

**Not here: the OCTA (Orange County) bus cycle.** It is the cycle I would have
reached for first, and the one the WVU transit-bus literature is built around,
but I could not find a freely redistributable second-by-second trace — DriveCAT
was unreachable and the rest were paywalled or summary-only. Rather than
digitise a figure from a paper and present it as OCTA data, the transit case is
carried by CBD as a close functional substitute, with Manhattan and NY Bus for
the severe end. If a proper trace turns up, adding a `case 'octa'` to
`cycle_mph10` is the entire integration effort.

**A caveat for the highway case.** UDDS and HWFET are *light-duty* chassis
cycles. Their accelerations are car accelerations, and a 17.7 t bus asked to
follow them demands more than its rated traction power on a handful of samples;
`ebus_range all` flags those with `[!]` rather than silently clipping the trace.
They are used for their speed profiles, not as certified bus cycles. The
genuinely coach-relevant highway result is the steady-cruise sweep in
`ebus_range sweeps`, a constant-speed operating point rather than a drive
cycle.

## Feeding bms-sim

The pack current this model produces is the stimulus my
[bms-sim](https://github.com/) project needs for its drive-cycle milestone —
two projects, one dataset.

```bash
matlab -batch "ebus_range profiles"
```

writes `./profiles/*.csv` as `time_s, current_A` (positive = discharge). The
conversion assumption is stated in `export_bms_profile` and repeated in every CSV
header: **pack voltage is held at its nominal 650 V**, so `I = P / V_nom`. A
real pack sags under load and rises under regen, so this understates current at
high discharge by roughly the voltage error — order 5–10 % at the extremes.
That is the right trade here, because bms-sim models the voltage response
itself; deriving current from a voltage it is about to compute would be
circular. This project owns the road load; that one owns the electrochemistry.

## Limitations

- **Quasi-static.** The speed trace is taken as given. No motor efficiency map
  (η_dt is one constant), no torque-speed envelope, no gear selection. A real
  motor is well below 90 % at low load.
- **It will not tell you a trace is impossible.** When a cycle demands more than
  rated traction power the model reports it (the ⚠ rows) but still follows the
  trace. UDDS and HWFET are light-duty chassis cycles — their accelerations are
  car accelerations. They are used here for their speed profiles, not as
  certified bus cycles.
- **No battery model.** Constant nominal voltage, no sag, no internal
  resistance, no temperature, no SOC-dependent power limits. That is bms-sim's
  job.
- **No thermal or HVAC model.** Auxiliaries are a constant kW.
- **Grade is a scalar or a hand-made vector**, not a real route elevation
  profile.
- **Regen recovery is probably optimistic.** η_regen = 0.85 assumes blended
  braking always prefers the motor. Real drivers and real brake blending are
  less disciplined; the `eta_regen` sweep exists to bound this.

## Possible extensions

- A real elevation profile for a specific route, instead of a scalar grade.
- A motor efficiency map replacing the single η_dt constant.
- Cold-weather operation: cabin heating raises the auxiliary load several-fold,
  and because auxiliaries are drawn against time, the cost per km is worst on
  exactly the slow stop-and-go routes transit buses run.
- Rebuilding it in Simulink — and in that order, so the validated MATLAB version
  is the reference to check the Simulink one against.

## Layout

Two files by design: this README and `ebus_range.m`. The single file holds the
model, the five published drive cycles, the plots and the test suite, so there
is no path setup and nothing to install.

The cost of that choice: the drive cycles are embedded as numeric literals
rather than CSVs you can open in a spreadsheet, the `data/raw` EPA and ADVISOR
sources are gone along with the converter that produced the CSVs from them, and
the generated figures are no longer committed — run `ebus_range all` to get
them back.

## License

MIT.

```
MIT License

Copyright (c) 2026 Karim

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
