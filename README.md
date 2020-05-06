# Fitting-density-profiles-for-disk-galaxies
This code performs different analysis on isolated disk galaxies:
0. Read gadget/gizmo output file into ascii
1. Particle selection and fitting of surface and vertical density profiles: for spherical components (bulge and halo) and for disky components.
2. Estimation of Star formation rates (SFR) based on either gas particles or stellar particles
3. Magnetic field stimation in case of Magneto-Hydrodynamic simulations.
4. Estimation of stability quatities such us dispersion velocities and rotational velocities.

Using the code: See param file
Compile the code: See Makefile options
runing the code: ./disk_surface_density.x SNAP param-file
