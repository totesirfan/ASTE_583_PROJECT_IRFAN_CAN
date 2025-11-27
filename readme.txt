ASTE 583: Lunar Trailblazer Navigation Project
Contributors: Irfan and Can
Date: Fall 2025

========================================================================
1. PROJECT OVERVIEW
========================================================================
This codebase implements orbit determination and maneuver planning for 
the Lunar Trailblazer mission recovery scenario. It includes high-fidelity
heliocentric dynamics (Sun + Earth + SRP), radiometric measurement 
modeling, and trajectory analysis.

Current Status: Phase 1 (Reference Trajectory & Model Validation) Complete.

========================================================================
2. LOCAL SETUP INSTRUCTIONS (READ FIRST)
========================================================================
To run this code on a new machine (e.g., Can's computer), follow these steps:

1. PREREQUISITES:
   - MATLAB (R2023b or newer recommended)
   - SPICE Toolkit for MATLAB (MICE) downloaded from NAIF.
   - SPICE Kernels downloaded from NAIF:
     * de440s.bsp (Planetary Ephemeris)
     * naif0012.tls.pc (Leap Seconds)

2. CONFIGURATION:
   - Open 'init_project.m'.
   - UPDATE line 6: Set 'mice_root' to the path where you extracted MICE.
   - UPDATE lines 16-19: Set the paths to where you saved the .bsp and .tls kernels.

3. INITIALIZATION:
   - Run 'init_project' in the Command Window before running any other script.
   - If successful, it will print: "SUCCESS: MICE added to MATLAB path."

========================================================================
3. FILE STRUCTURE
========================================================================
The codebase is organized into "Libraries" (Core Logic) and "Phases" 
(Execution Scripts) to separate physics from analysis.

[LIBRARIES] - DO NOT MODIFY UNLESS NECESSARY
- lib_constants.m       : Single source of truth for physics constants, 
                          station coordinates, and initial states.
- lib_dynamics.m        : Equations of Motion (Sun Gravity + Earth 3rd 
                          Body + SRP).
- lib_measurements.m    : Computes Range/Doppler observations and the 
                          H-matrix (Partials) for estimation.

[PHASE 1: MODELING & REFERENCE TRAJECTORY]
- Phase1_Test_Integrity.m : Validates the math by comparing analytical 
                            partials (H-matrix) against numerical finite 
                            differences. (Must PASS before proceeding).
- Phase1_Main_RefTraj.m   : Propagates the nominal trajectory, applies 
                            the LTM maneuver, and plots the Lunar Flyby
                            in Heliocentric and Geocentric views.
- Phase1_Viz_GroundTrack.m: Generates the reference ground track map.

[PHASE 2: ORBIT DETERMINATION] (Upcoming)
- Phase2_Load_Data.m      : Parses CSV measurement files.
- Phase2_Batch_Filter.m   : Least Squares estimator.

========================================================================
4. EXECUTION ORDER
========================================================================
1. >> init_project
   (Loads SPICE. Do this once per MATLAB session).

2. >> Phase1_Test_Integrity
   (Verifies that physics and math are correct. Max Error should be < 1e-4).

3. >> Phase1_Main_RefTraj
   (Generates the Reference Trajectory and Flyby plots).

4. >> Phase1_Viz_GroundTrack
   (Generates the Ground Track map).

========================================================================
5. NOTES FOR COLLABORATORS
========================================================================
- Always run 'Phase1_Test_Integrity' after modifying any Library file.
- 'lib_measurements.m' contains the logic for the "Bias Switch" (Day 6).
- All coordinate frames are standardized to Sun-Centered EMO2000.
