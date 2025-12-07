# ASTE 583 Project – Lunar Trailblazer Navigation

**Repository:** `totesirfan-aste_583_project_irfan_can`  
**Course:** ASTE 583 – Space Navigation Principles and Practice
**Group Members:**
-Irfan 
-Can

This repository implements the full Orbit Determination (OD) and maneuver analysis pipeline for the **Lunar Trailblazer** class project. It adheres to the **official ASTE 583 lecture slides** and the **ASTE583_Project.pdf handout**.

## Project Overview

The workflow includes:
- **Phase 1:** Prime-nav reference trajectory and ground track generation.
- **Phase 2 (Prelim):** Preliminary OD (0–6 days) using Batch Least Squares and Sequential Kalman Filter (CKF).
- **Phase 2 (Final):** Final OD (0–14 days) using an iterated CKF at the detection epoch.
- **Phase 3:** Cleanup maneuver statistics (LTM + LCM Monte Carlo) and LTM maneuver execution signature analysis.

All dynamics and measurement models are consistent with the project handout and the ASTE 583 slide decks.

---

## 1. Repository Structure

```text
totesirfan-aste_583_project_irfan_can/
├── init_project.m                     # Path setup & kernel loading
├── lib_constants.m                    # Centralized constants & configuration
├── lib_dynamics.m                     # Dynamics model (Sun, Earth, SRP)
├── lib_measurements.m                 # Measurement model (Range, Doppler)
├── load_project_measurements.m        # Data loader for CSVs
├── Phase1_Main_RefTraj.m              # Phase 1: Reference Trajectory
├── Phase2_Prelim_Batch.m              # Phase 2: Prelim Batch LS
├── Phase2_Prelim_Kalman.m             # Phase 2: Prelim CKF
├── Phase2_Prelim_Plots.m              # Phase 2: Prelim Residual Plots
├── Phase2_Prelim_Passthru.m           # Phase 2: Passthru test (6-14d)
├── Phase2_Final_Kalman.m              # Phase 2: Final Iterated CKF
├── Phase2_FinalKalman_PostProcess.m   # Phase 2: Final Residuals & RTN
├── Phase3_CleanupStats.m              # Phase 3: Monte Carlo Statistics
└── Phase3_ManeuverExecution.m         # Phase 3: Maneuver RR Signature
```

---

## 2. External Dependencies

### MATLAB
- Tested in MATLAB (R202x).
- **Core functions:** `ode45`, `readtable`, `interp1`, basic linear algebra (`chol`, `eig`, etc.).

### SPICE / MICE
This project uses NAIF SPICE via MICE for ephemerides and time conversions. You must edit `init_project.m` to point to your local installation:

```matlab
mice_root = 'C:\Users\irfan\OneDrive\Documents\MATLAB\mice';
kernels = {
    'C:\Users\irfan\OneDrive\Documents\MATLAB\kernels\de440s.bsp';
    'C:\Users\irfan\OneDrive\Documents\MATLAB\kernels\naif0012.tls.pc'
};
```
**Required Kernels:**
- `de440s.bsp`: Planetary ephemerides
- `naif0012.tls.pc`: Leapseconds / time conversion

### Measurement Data
The OD pipeline expects two CSV files in the MATLAB path:
1. `ASTE583_Project_LTB_Measurements_0-6D_Truth.csv`
2. `ASTE583_Project_LTB_Measurements_6D-14D_Truth.csv`

---

## 3. Frames, Units & Conventions

- **Frame:** Sun-centered EMO2000 (Earth Mean Obliquity of J2000).
- **Earth Rotation:** Simple rotation about EMO z-axis using Greenwich sidereal angle at detection epoch (`const.phi_G_detect`).
- **State Vector (10-state OD):**
  ```matlab
  X = [r_x, r_y, r_z, v_x, v_y, v_z, k_SRP, rho_bias, lat_4, lon_4]'
  ```
- **Units:**
  - Position: km
  - Velocity: km/s
  - Time: Seconds (ET) and days since detection
  - Station Lat/Lon: Radians
  - Doppler bias: km/s

---

## 4. Quick Start

### 4.1 One-Time Setup
Run this once per MATLAB session to load MICE:
```matlab
init_project();
const = lib_constants();
```

### 4.2 Typical Workflow (Recommended Order)

**Phase 1: Reference Trajectory**
```matlab
Phase1_Main_RefTraj;
% Output: Trajectory plots, ground track, closest approach stats.
```

**Phase 2: Preliminary OD (0–6 days)**
```matlab
Phase2_Prelim_Batch;   % Batch LS -> ASTE583_PrelimBatch_Results.mat
Phase2_Prelim_Kalman;  % CKF -> ASTE583_PrelimKalman_Results.mat
Phase2_Prelim_Plots;   % Compare residuals and 3-sigma ellipses
Phase2_Prelim_Passthru;% Check generalization to 6-14d data
```

**Phase 2: Final OD (0–14 days)**
```matlab
% Recommended: use 'handout' prior for physical consistency
Phase2_Final_Kalman('handout'); 
Phase2_FinalKalman_PostProcess; 
% Output: ASTE583_FinalKalman_Results.mat, residuals, RTN uncertainties.
```

**Phase 3: Maneuvers & Statistics**
```matlab
Phase2_Cleanup_Statistics;            % Monte Carlo cleanup stats (defined in Phase3_CleanupStats.m)
Phase3_LTM_ManeuverExecution_Final;   % LTM RR signature (defined in Phase3_ManeuverExecution.m)
```

---

## 5. Script-by-Script Overview

### 5.1 Core Utilities
- **`init_project.m`**: Adds MICE to path and furnishes kernels.
- **`lib_constants.m`**: Returns `const` struct with body parameters, station coords, time constants, and maneuver definitions.
- **`lib_dynamics.m`**: Dynamics model including Sun point-mass, Earth third-body, and Cannonball SRP. Supports State+STM propagation.
- **`lib_measurements.m`**: Measurement model (Y) and H-matrix. Includes time-limited Doppler bias (active 0–6 days only).
- **`load_project_measurements.m`**: Merges and sorts the two source CSV files.

### 5.2 Phase 1
- **`Phase1_Main_RefTraj.m`**: Propagates prime-nav trajectory. Generates heliocentric plots, geocentric flyby views, and global ground tracks with station markers.

### 5.3 Phase 2 (Preliminary)
- **`Phase2_Prelim_Batch.m`**: 10-state Batch LS over 0–6 days. Solves for state, k_SRP, bias, and Station 4 coordinates.
- **`Phase2_Prelim_Kalman.m`**: Single-pass CKF over 0–6 days.
- **`Phase2_Prelim_Plots.m`**: Compares Batch vs. Kalman residuals and plots 3σ RTN ellipses at detection epoch.
- **`Phase2_Prelim_Passthru.m`**: Propagates prelim solutions forward to check fit against 6–14 day data.

### 5.4 Phase 2 (Final)
- **`Phase2_Final_Kalman.m`**: Iterated CKF at detection epoch using full 0–14 day arc. Supports 'handout', 'batch', or 'kalman' priors.
- **`Phase2_FinalKalman_PostProcess.m`**: Generates final residual plots (mm/s) and propagates covariance to LTM to compute RTN uncertainties.

### 5.5 Phase 3
- **`Phase3_CleanupStats.m`**: Defines `Phase2_Cleanup_Statistics()`. Performs Monte Carlo analysis (default N=10,000) to estimate cleanup ΔV (LCM) and total ΔV budget compliance.
- **`Phase3_ManeuverExecution.m`**: Defines `Phase3_LTM_ManeuverExecution_Final()`. Visualizes the LTM maneuver signature in range-rate for all ground stations.

---

## 6. Implementation Notes
- **Function Names:** - `Phase3_CleanupStats.m` defines the function `Phase2_Cleanup_Statistics()`.
    - `Phase3_ManeuverExecution.m` defines the function `Phase3_LTM_ManeuverExecution_Final()`.
- **Integrator:** `ode45` is used with tight tolerances (`RelTol ~ 1e-10`, `AbsTol ~ 1e-9`) for numerical stability.
- **RNG:** Monte Carlo simulations use `rng(583)` for reproducibility.

---

## 7. Credits
**Authors:** Irfan & Can  
**Course:** ASTE 583 – Navigation & Estimation (USC)  
**Reference:** Based on the Lunar Trailblazer OD project handout and ASTE 583 lecture slides.
