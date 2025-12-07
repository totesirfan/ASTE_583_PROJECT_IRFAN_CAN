# ASTE 583 Project – Lunar Trailblazer Navigation

**Repository:** `totesirfan-aste_583_project_irfan_can`
**Course:** ASTE 583 – Space Navigation Principles and Practice
**Group Members:** Irfan & Can

This repository implements the full Orbit Determination (OD) and maneuver analysis pipeline for the **Lunar Trailblazer** class project. It adheres strictly to the **official ASTE 583 lecture slides** and the **ASTE583_Project.pdf** handout.

---

## Project Overview

The workflow includes:

- **Phase 1:** Prime-nav reference trajectory and ground track generation.
- **Phase 2 (Prelim):** Preliminary OD (0–6 days) using Batch Least Squares and a Sequential Kalman Filter (CKF).
- **Phase 2 (Final):** Final OD (0–14 days) using an iterated CKF at the detection epoch.
- **Phase 3:** Cleanup maneuver statistics (LTM + LCM Monte Carlo) and LTM maneuver execution signature analysis.

All dynamics and measurement models are consistent with the project handout and the ASTE 583 slide decks.

---

## 1. Repository Structure

```text
totesirfan-aste_583_project_irfan_can/
├── init_project.m                     # Path setup & SPICE kernel loading
├── lib_constants.m                    # Centralized constants & configuration
├── lib_dynamics.m                     # Dynamics model (Sun, Earth, SRP)
├── lib_measurements.m                 # Measurement model (Range, Doppler)
├── load_project_measurements.m        # Data loader for CSVs
├── Phase1_Main_RefTraj.m              # Phase 1: Reference Trajectory
├── Phase2_Prelim_Batch.m              # Phase 2: Prelim Batch LS
├── Phase2_Prelim_Kalman.m             # Phase 2: Prelim CKF
├── Phase2_Prelim_Plots.m              # Phase 2: Prelim residual plots & RTN
├── Phase2_Prelim_Passthru.m           # Phase 2: Passthru test (6–14 d)
├── Phase2_Final_Kalman.m              # Phase 2: Final Iterated CKF (0–14 d)
├── Phase2_FinalKalman_PostProcess.m   # Phase 2: Final residuals & LTM RTN
├── Phase3_ManeuverExecution.m         # Phase 3: LTM maneuver RR signature
└── Phase3_CleanupStats.m              # Phase 3: Cleanup Monte Carlo statistics
```

---

## 2. External Dependencies

### MATLAB

* Tested in MATLAB R202x.
* Uses only base MATLAB functionality:
  * ODE integration: `ode45`
  * Table I/O: `readtable`
  * Interpolation: `interp1`
  * Linear algebra: `chol`, `eig`, etc.

No toolboxes beyond base MATLAB are required.

### SPICE / MICE

This project uses **NAIF SPICE via MICE** for ephemerides and time conversions.
You **must** edit `init_project.m` to point to your local installation:

```matlab
% In init_project.m
mice_root = 'C:\path\to\mice';
kernels = {
    'C:\path\to\kernels\de440s.bsp';
    'C:\path\to\kernels\naif0012.tls.pc'
};
```

**Main MICE routines used:**

* `cspice_str2et` – convert UTC strings to ET
* `cspice_spkezr` – state vectors from SPK kernels
* `cspice_et2utc` – convert ET back to human-readable UTC strings

If any of these are “undefined”, your MICE path is not set correctly.

### Measurement Data

The OD pipeline expects two CSV files either in the **current MATLAB folder** or **on the MATLAB path**:

1. `ASTE583_Project_LTB_Measurements_0-6D_Truth.csv`
2. `ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv`

`load_project_measurements.m` loads, merges, and time-sorts these into a single structure used by all Phase 2 scripts.

---

## 3. Frames, Units & Conventions

* **Frame:** Sun-centered **EMO2000** (Earth Mean Obliquity of J2000).
* **Earth Rotation:** Simple rotation about the EMO z-axis using the Greenwich sidereal angle at detection epoch (`const.phi_G_detect`) and a constant rotation rate `const.we`.
* **State Vector (10-state OD):**

  ```matlab
  X = [r_x, r_y, r_z, v_x, v_y, v_z, k_SRP, rho_bias, lat_4, lon_4]';
  ```
* **Units:**
  * Position: km
  * Velocity: km/s
  * Time: ET seconds internally, days since detection in plots
  * Station latitude / longitude: radians
  * Doppler bias: km/s

All modeling choices (SRP cannonball, third-body Earth, use of EMO2000, etc.) follow the ASTE 583 slides and project handout.

---

## 4. Quick Start

### 4.1 One-Time Setup (per MATLAB session)

1. Set your working directory to the repo folder.
2. Run:

```matlab
init_project();       % Adds MICE to path and furnishes kernels
const = lib_constants();  % Optional: inspect constants & epochs
```

If this fails with “undefined function `cspice_spkezr`” or missing kernel warnings, fix the paths in `init_project.m` before proceeding.

### 4.2 Recommended Execution Order

#### Phase 1: Reference Trajectory

```matlab
Phase1_Main_RefTraj;
% Output:
%  - Heliocentric trajectory plot
%  - Geocentric + lunar flyby view
%  - Global ground track with stations
%  - Printed closest approach distance, altitude, and UTC
```

#### Phase 2: Preliminary OD (0–6 days)

```matlab
Phase2_Prelim_Batch;    % Batch LS -> ASTE583_PrelimBatch_Results.mat
Phase2_Prelim_Kalman;   % CKF      -> ASTE583_PrelimKalman_Results.mat
Phase2_Prelim_Plots;    % Residuals & RTN 3σ ellipses at t0
Phase2_Prelim_Passthru; % Check generalization on 6–14 d data
```

#### Phase 2: Final OD (0–14 days)

```matlab
% Recommended: use 'handout' prior for a physically consistent solution
Phase2_Final_Kalman('handout');
Phase2_FinalKalman_PostProcess;
% Output:
%  - ASTE583_FinalKalman_Results.mat
%  - Final prefit/postfit residual plots (mm/s, km)
%  - LTM RTN covariance and 1σ / 3σ position uncertainties
```

> **Note:** The code also supports `'batch'` and `'kalman'` as priors for `Phase2_Final_Kalman`, but these can drive the SRP scale factor `k_SRP` to an unphysical negative value over the full 0–14 day arc. For the final reported solution, we use the **'handout'** prior and stop the CKF iterations before entering that unphysical basin.

#### Phase 3: Maneuvers & Statistics

```matlab
Phase3_ManeuverExecution;   % LTM range-rate signature for all stations
Phase3_CleanupStats;        % Monte Carlo cleanup statistics (LTM + LCM)
```

Make sure the function names you call here match the primary function names in `Phase3_ManeuverExecution.m` and `Phase3_CleanupStats.m` (MATLAB requires the main function name to match the filename).

---

## 5. Outputs & Artifacts

| Script / Function                | Key Outputs                                                               |
| -------------------------------- | ------------------------------------------------------------------------- |
| `Phase1_Main_RefTraj`            | Figures (Sun/Earth/Moon views, ground track) + printed CA distance & UTC  |
| `Phase2_Prelim_Batch`            | `ASTE583_PrelimBatch_Results.mat` (X₀, P₀, per-iteration diagnostics)      |
| `Phase2_Prelim_Kalman`           | `ASTE583_PrelimKalman_Results.mat` (X₀, P₀, pre/post RR residuals)         |
| `Phase2_Prelim_Plots`            | Plots of Batch vs Kalman RR residuals, RTN 3σ ellipses at detection epoch |
| `Phase2_Prelim_Passthru`         | Plots of 6–14 d range / Doppler passthru residuals (Batch vs Kalman)      |
| `Phase2_Final_Kalman`            | `ASTE583_FinalKalman_Results.mat` (final 10-state X₀, P₀, full residuals) |
| `Phase2_FinalKalman_PostProcess` | Final RR / range residual plots + printed LTM RTN σ_R, σ_T, σ_N           |
| `Phase3_ManeuverExecution`       | Plots of RR_noLTM vs RR_LTM and ΔRR around LTM for all stations           |
| `Phase3_CleanupStats`            | Histograms of LCM dispersion & ΔV, printed ΔV₉₉ and dV_budget comparison  |

These are the files and figures referenced in the project write-up.

---

## 6. Script-by-Script Overview

### 6.1 Core Utilities

* **`init_project.m`**
  Adds MICE to the MATLAB path and furnishes the required SPICE kernels. Must be run at the start of each MATLAB session.

* **`lib_constants.m`**
  Returns a struct `const` with:
  * Gravitational parameters, radii, SRP constants, and EMO rotation matrix.
  * Station coordinates (lat/lon), detection epoch, GST at detection.
  * Maneuver definitions (LTM and LCM epochs, nominal ΔV, execution σ).
  * A priori covariance matrices for 6-state and 10-state OD.

* **`lib_dynamics.m`**
  Implements the Sun-centered EMO2000 dynamics:
  * Sun point-mass gravity,
  * Earth third-body gravity,
  * Cannonball SRP scaled by `k_SRP`.
    Supports both:
  * State-only propagation, and
  * Augmented state + STM propagation (10-state + 10×10 STM).

* **`lib_measurements.m`**
  Implements the radiometric measurement model:
  * Two-way range and range-rate (Doppler),
  * Time-limited Doppler bias active only over the first 6 days,
  * Optional estimation of station 4 latitude/longitude,
  * Returns both the computed measurement `Y_comp` and the linearized H-matrix.

* **`load_project_measurements.m`**
  Loads the two CSV files, merges them into a single time-ordered data set, and returns a struct with:
  * `time_sec`, `station_id`, `range_km`, `rr_kmps`.

### 6.2 Phase 1

* **`Phase1_Main_RefTraj.m`**
  Propagates the prime-nav reference trajectory from detection to 251 days, including the LTM maneuver. Produces:
  * Heliocentric trajectory plots,
  * Geocentric + lunar flyby plots,
  * Global ground track with station markers,
  * Closest-approach distance, altitude, and UTC time.

User-adjustable time steps (`dt_state_sec`, `dt_track_sec`) are defined at the top of the script.

### 6.3 Phase 2 – Preliminary OD (0–6 days)

* **`Phase2_Prelim_Batch.m`**
  Performs 10-state Batch Least Squares OD over the 0–6 day arc. Solves for:
  * Initial state [r₀, v₀],
  * SRP scale `k_SRP`,
  * Doppler bias,
  * Antarctica station latitude/longitude.
    Saves all results and per-iteration residuals to `ASTE583_PrelimBatch_Results.mat`.

* **`Phase2_Prelim_Kalman.m`**
  Runs a single-pass CKF over 0–6 days using the same 10-state model. Produces a posterior state/covariance at detection epoch (`X0_KF`, `P0_KF`) and Doppler prefit/postfit residual histories.

* **`Phase2_Prelim_Plots.m`**
  Uses the batch and Kalman results to:
  * Plot RR prefit/postfit residuals (in mm/s), and
  * Compare 3σ RTN ellipses at detection epoch.

* **`Phase2_Prelim_Passthru.m`**
  Propagates both prelim solutions forward and evaluates residuals on the 6–14 day measurements to assess generalization.

### 6.4 Phase 2 – Final OD (0–14 days)

* **`Phase2_Final_Kalman.m`**
  Performs an **iterated CKF at detection epoch** using the full 0–14 day radiometric data set.
  Supports three choices of prior:
  * `'handout'` – a priori from the project handout (recommended for final solution),
  * `'batch'` – a priori from the prelim Batch OD,
  * `'kalman'` – a priori from the prelim Kalman OD.

  Internally:
  * Propagates the nominal 10-state and STM from detection to the last measurement,
  * Forms measurement sensitivities using `H * Φ`,
  * Applies sequential CKF updates to a single ΔX₀,
  * Iterates until corrections are small or the maximum iteration count is reached.

* **`Phase2_FinalKalman_PostProcess.m`**
  Post-processes the final CKF result by:
  * Plotting prefit/postfit Doppler and range residuals (per station, in mm/s and km),
  * Propagating the 10×10 covariance to the LTM epoch,
  * Transforming into the RTN frame and printing 1σ / 3σ values for the radial, along-track, and cross-track position components.

These outputs are used directly to answer the final OD and LTM targeting questions in the handout.

### 6.5 Phase 3 – Maneuvers & Statistics

* **`Phase3_ManeuverExecution.m`**
  Computes the predicted range-rate signature of the LTM maneuver:
  * Propagates both a “no-LTM” and “with-LTM” trajectory around the LTM epoch,
  * Evaluates range-rate for all four stations,
  * Plots RR_noLTM vs RR_LTM and ΔRR (LTM − noLTM) in mm/s.

* **`Phase3_CleanupStats.m`**
  Performs the cleanup-maneuver Monte Carlo analysis using the **final CKF a posteriori covariance**:
  * Draws random 10-state samples at detection epoch,
  * Propagates through LTM with execution error and to LCM,
  * Computes position/velocity dispersion at LCM,
  * Applies a linear feedback law at LCM to target the 30-day-post-LCM reference,
  * Records LCM ΔV, LTM ΔV, and total ΔV for each trial,
  * Outputs ΔV statistics including ΔV₉₉ and the fraction of trials within the project ΔV budget.

---

## 7. Implementation Notes & Troubleshooting

* **Integrator Settings:**
  All ODEs are integrated with tight tolerances (`RelTol ≈ 1e-10`, `AbsTol ≈ 1e-9`) for numerical stability over the long arcs.

* **Random Seeds:**
  Monte Carlo runs use `rng(583)` for reproducibility.

* **SRP Scale Factor Behavior:**
  When using `'batch'` or `'kalman'` priors (or a forward EKF over 0–14 days), the unconstrained least-squares solution tends to drive `k_SRP` to a negative value. This is mathematically optimal for the chosen force model but **physically impossible**; the final reported solution therefore uses the **'handout'** prior and an iterated CKF that stays in the physically meaningful basin (`k_SRP > 0`).

* **Common Errors:**

  * `Undefined function or variable 'cspice_spkezr'`
    → MICE is not on the MATLAB path. Edit `init_project.m` and re-run `init_project()`.
  * `Unable to open ASTE583_Project_*.csv`
    → Ensure the measurement CSVs are in the current directory or on the MATLAB path.
  * Strange ground track or station geometry
    → Check `const.phi_G_detect`, station lat/lon, and that MICE kernels are being furnished correctly.

---

## 8. Credits

**Authors:** Irfan & Can
**Course:** ASTE 583 – Navigation & Estimation (USC)

This implementation follows the Lunar Trailblazer OD project handout and the ASTE 583 lecture slides as the authoritative references.
