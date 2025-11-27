ASTE 583: Lunar Trailblazer Navigation Project
Contributors: Irfan and Can
Date: Fall 2025

========================================================================
1. PROJECT OVERVIEW
========================================================================
This codebase implements orbit determination and maneuver planning for 
the Lunar Trailblazer mission recovery scenario. It includes high-fidelity
heliocentric dynamics, radiometric measurement modeling, and statistical
maneuver analysis.

========================================================================
2. SETUP INSTRUCTIONS
========================================================================
1. Ensure the 'mice' toolkit folder is present.
2. Ensure the 'kernels' folder contains 'de440s.bsp' and 'naif0012.tls.pc'.
3. Open 'init_project.m' and update the paths if necessary.
4. Run 'init_project.m' to load SPICE kernels into MATLAB memory.

========================================================================
3. FILE STRUCTURE
========================================================================
The code is organized into Libraries (Core Logic) and Phases (Execution).

[LIBRARIES]
- lib_constants.m       : Single source of truth for physics constants, 
                          station coordinates, and initial states.
- lib_dynamics.m        : Equations of Motion (Sun Gravity + Earth 3rd 
                          Body + SRP).
- lib_measurements.m    : Computes Range/Doppler observations and the 
                          H-matrix (Partials) for estimation.

[PHASE 1: MODELING & REFERENCE]
- Phase1_Test_Integrity.m : Validates the analytical H-matrix against 
                            numerical finite differences. (MUST PASS)
- Phase1_Main_RefTraj.m   : Propagates the nominal trajectory, applies 
                            the LTM maneuver, and plots the Lunar Flyby.
- Phase1_Viz_GroundTrack.m: Generates the "Dot Cloud" ground track map 
                            for the reference trajectory.

[PHASE 2: ORBIT DETERMINATION] (Upcoming)
- Phase2_Load_Data.m      : Parses CSV measurement files.
- Phase2_Batch_Filter.m   : Least Squares estimator.

========================================================================
4. EXECUTION ORDER
========================================================================
1. Run 'init_project' (Once per session).
2. Run 'Phase1_Test_Integrity' to verify physics models.
3. Run 'Phase1_Main_RefTraj' to generate trajectory plots.
4. Run 'Phase1_Viz_GroundTrack' to generate the ground track figure.

========================================================================
