## Introduction to ROMS

### Background
The Regional Ocean Modeling System (ROMS) is a well-designed Fortran package, which simulates the free-surface geophysical fluid dynamics system, using an hydrostatic, primitive equatesiton with Boussinesq approximation. ROMS is developed with the Nonlinear integration kernel which has a wide application in nonlinear fluid dynamic studies.terrain-following vertical coordinate is applied in ROMS, in order to achieve better vertical resolution in shallow water and areas with complex bathymetry.

###Dynamic Equations

The governing dynamical equations of three-dimensional, free-surface, Reynolds-averaged Navier-Stokes equations are

$$
\frac{\partial H_zu}{\partial t} + \frac{ \partial (uH_zu)}{\partial x} + \frac{ \partial( v H_z u)}{\partial y}  + \frac{\partial(\Omega H_zu)}{\partial s}  - f H_z v = - \frac{ H_z}{\rho_0}\frac{ \partial p}{\partial x} -   H_z g\frac{ \partial \zeta }{\partial x} - \frac{\partial}{\partial s}(-K_M \frac{\partial u}{\partial z} - \frac{v}{H_z}\frac{\partial u}{\partial s})
$$

$$
\frac{\partial H_zv}{\partial t} + \frac{ \partial (uH_zv)}{\partial x} + \frac{ \partial( v H_z v)}{\partial y}  + \frac{\partial(\Omega H_zv)}{\partial s}  + f H_z u = - \frac{ H_z}{\rho_0}\frac{ \partial p}{\partial y} -   H_z g\frac{ \partial \zeta }{\partial y} - \frac{\partial}{\partial s}( -K_M \frac{\partial v}{\partial z} - \frac{v}{H_z}\frac{\partial v}{\partial s})
$$

with the continuity equation

$$
\frac{\partial \zeta}{\partial t} + \frac{ \partial (H_zu)}{\partial x} + \frac{ \partial( H_z v)}{\partial y}  + \frac{\partial(H_z \Omega)}{\partial s}  = 0
$$

equation of State

$$
\rho = \rho (T, S, P)
$$
and scalar transport equation for temperature and salinity

$$
\frac{\partial (H_z C)}{\partial t} + \frac{ \partial (u H_z C)}{\partial x} + \frac{ \partial( v H_z C)}{\partial y}  + \frac{\partial(\Omega H_z C)}{\partial s} = -\frac{\partial}{\partial s}( -K_C \frac{\partial C}{\partial z} - \frac{\upsilon_\theta}{H_z}\frac{\partial C}{\partial s})
$$

Here, all variables are:


| name                                | Description                             |
| ----------------------------------- | --------------------------------------- |
| $u$                                 | horizontal velocity in x direction      |
| $v$                                 | horizontal velocity in y direction      |
| $s$                                 | the scaled sigma coordinate             |
| $\Omega$                            | vertical velocity in (sigma coordinate) |
| $\zeta$                             | free-surface elevation                  |  

| name                                | Description                             |
| ----------------------------------- | --------------------------------------- |
| $H_z$                               | vertical stretching factor              |
| $f$                                 | Coriolis parameter                      |
| $K_M$, $K_C$                        | vertical eddy viscosity and diffusivity |
| $p$                                 | pressure                                |
| $\rho_0$                            | reference density                       |
| $g$                                 | acceleration due to gravity             |
| $C$                                 | tracer (Temperature, Salinity, etc.)    |
| Spatial Discretization (horizontal) | Square Grid                             |
| Spatial Discretization (vertical)   | Terrain Coordinate                      |
| Time Discretization                 | Time Step                               |

In the Navier-Stokes equations, the hydrostatic approximation is used since the horizontal scale is usually larger than the vertical scale thus making the approximation valid. The Boussinesq approximation is also applied here, ignoring the density difference except the terms arisen from the gravity.


### Vertical S-coordinate

In ROMS, the terrain following coordinate system is used in vertical direction, which means that the number of layers for the ocean is the same everywhere while the thickness of each layer varies with the bathymetry of the specific location. There are several options for the transformation equations built in ROMS, controlled by two input argument, Vtransform and Vstretch. For our project, we choose Vtransform = 2 and Vstretch = 4. The details of these two choices will be discussed below.
#### Transformation Equations
Vtransform being 2 means that we are using the formulation developed by A. Shchepetkin in 2005 [1]:

$$
z(x,y,\sigma,t) = \zeta(x,y,t)+[\zeta(x,y,t)+h(x,y)]S(x,y,\sigma)
$$

$$
S(x,y,\sigma) = \frac{h_c\sigma+h(x,y)C(\sigma)}{h_c+h(x,y)}
$$
where $h_c$ is a critical depth controlling the resolution and stretching, which will be desctribed in detail below, $\zeta$ is the free surface elevation, $\sigma$ is a fractional stretching coordinate from $-1\leq\sigma\leq0$ (-1 for the ocean bottom and 0 for the sea surface), and C($\sigma$) is a nondimensinal stretching function ranging from $-1\leq C(\sigma)\leq0$.
It is convenient to define the vertical stretching factor as
$$
H_z\equiv\frac{\partial z}{\partial\sigma}
$$
Then $H_z(x,y,\sigma,t)$ is the vertical grid thickness. In ROMS, it is computed as $\Delta z/\Delta\sigma$.
####Stretching Functions
(Ref. ROMS Manual)For the option Vstretch = 4, the stretching function is defined as a double stretching function:
Surface refinement function as
$$
C(\sigma) = \frac{1-cosh(\theta_s\sigma)}{cosh(\theta_s)-1}\qquad      for\ \theta_s> 0\\

C(\sigma) = -\sigma^2\qquad      for\ \theta_s\leq0
$$
Bottom refinement function as
$$
C(\sigma) = \frac{e^{\theta_BC(\sigma)}-1}{1-e^{-\theta_B}}\qquad      for\ \theta_B>0
$$
The rage of the parameters are $0\leq\theta_s\leq 10$ and $0\leq\theta_B\leq 4$.

###Turbulence Closure
When the Reyhnolds-averaged Navier-Stokes equations are first derived, there are terms including the average in the perturbations of velocity, such as $\overline{u'w'}$, making the number of unknown variables greater than the number of equations. In order to solve the problem, the turbulence closure technique is used, and the equations are:
$$
\overline{u'w'} = -K_M\frac{\partial u}{\partial z}\\
\overline{v'w'} = -K_M\frac{\partial v}{\partial z}\\
\overline{C'w'} = -K_C\frac{\partial C}{\partial z}
$$

As shown above, these equations approximate those terms with the gradient of velocity/tracer, implying that they flow down the local gradient of u, v, C respectively. This method of turbulence closure is also called the Gradient Transport Theory or K-theory.

##Data Assimilation with ROMS

###Why Data Assimilation is Needed?

Study of the data assimilation with partial observation is necessary and a fundamental challenge meteorology and oceanography because, in practice, it is impossible to measure the exact variable states of the entire system, due to but not limited to the following reasons:

1. Making measurements will cost too much. For example, measuring the velocity fields in the deep ocean.

2. The measurements have errors because of the limits of equipments.

3. The model may be inaccurate, such as making inappropriate assumptions.

There are several different data assimilation method could be applied in the ROMS. Here the two commonly used ones are simple nudging method and incremental 4D-VAR (I4D-VAR).

###Simple Nudging

TBD

###I4D-VAR

TBD

##Appendix

###Details of the I4D-VAR Code (Reference or starting point for us to apply our dynamical nudging or DSPE method to the ROMS)

A module file mod_fourdvar.F is used to run the minimization of the cost function using Lanczos algorithm or descent algorithm. (More details are being studied in order for us to put in our dynamical nudging methos)

The cost function and its gradient is called to compute in the following section of code in is4dvar_ocean.h, which is also the part we are planning to modify to apply dynamical nudging method. The actual computation is done in tl_main3d.F.

```
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time-step tangent linear model: compute cost function.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  If first pass inner=0, initialize tangent linear state (increments,
!  deltaX) from rest. Otherwise, use trial initial conditions estimated
!  by the conjugate gradient algorithm in previous inner loop. The TLM
!  initial conditions are read from ITL(ng)%name, record 1.
!
          DO ng=1,Ngrids
            ITL(ng)%Rindex=1
!$OMP PARALLEL
            CALL tl_initial (ng)
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!  On first pass, initialize records 2, 3 and 4 of the ITL file to zero.
!
          IF (inner.eq.0.and.outer.eq.1) THEN
            DO ng=1,Ngrids
              CALL tl_wrt_ini (ng, LTLM1, Rec2)
              IF (exit_flag.ne.NoError) RETURN
              CALL tl_wrt_ini (ng, LTLM1, Rec3)
              IF (exit_flag.ne.NoError) RETURN
              CALL tl_wrt_ini (ng, LTLM1, Rec4)
              IF (exit_flag.ne.NoError) RETURN
            END DO
          END IF

#ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, activate writing and determine
!  output file name. The multiple file option is use to perturb initial
!  state and create ensembles.  The TLM final trajectory is written for
!  each inner loop on separated NetCDF files.
!
          DO ng=1,Ngrids
            LdefTLM(ng)=.TRUE.
            LwrtTLM(ng)=.TRUE.
            WRITE (TLM(ng)%name,10) TRIM(TLM(ng)%base), Nrun
          END DO
#endif
!
!  Activate switch to write out initial and final misfit between
!  model and observations.
!
          DO ng=1,Ngrids
            wrtMisfit(ng)=.FALSE.
            IF (((outer.eq.1).and.(inner.eq.0)).or.                     &
     &          ((outer.eq.Nouter).and.(inner.eq.Ninner))) THEN
              wrtMisfit(ng)=.TRUE.
            END IF
          END DO
!
!  Run tangent linear model. Compute misfit observation cost function,
!  Jo.
!
          DO ng=1,Ngrids
            IF (Master) THEN
              WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
            END IF
          END DO

!$OMP PARALLEL
#ifdef SOLVE3D
          CALL tl_main3d (RunInterval)
#else
          CALL tl_main2d (RunInterval)
#endif
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN

#ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, close current NetCDF file.
!
          DO ng=1,Ngrids
            IF (TLM(ng)%ncid.ne.-1) THEN
              SourceFile='is4dvar_ocean.h, ROMS_run'

              CALL netcdf_close (ng, iTLM, TLM(ng)%ncid)
              IF (exit_flag.ne.NoError) RETURN
            END IF
          END DO
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time step adjoint model backwards: compute cost function gradient.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize the adjoint model always from rest.
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            CALL ad_initial (ng)
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!  Time-step adjoint model backwards. The adjoint model is forced with
!  the adjoint of the observation misfit (Jo) term.
!
          DO ng=1,Ngrids
            IF (Master) THEN
              WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
            END IF
          END DO

!$OMP PARALLEL
#ifdef SOLVE3D
          CALL ad_main3d (RunInterval)
#else
          CALL ad_main2d (RunInterval)
#endif
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN
!
!  Clear adjoint arrays.  Is it needed?
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL initialize_ocean (ng, tile, iADM)
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
              CALL initialize_forces (ng, tile, iADM)
#endif
            END DO
!$OMP END PARALLEL
          END DO
!
```