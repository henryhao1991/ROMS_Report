## Introduction to ROMS

### Background
The Regional Ocean Modeling System (ROMS) is a well-designed Fortran package, which simulates the free-surface geophysical fluid dynamics system, using an hydrostatic, primitive equatesiton with Boussinesq approximation. ROMS is developed with the Nonlinear integration kernel which has a wide application in nonlinear fluid dynamic studies.terrain-following vertical coordinate is applied in ROMS, in order to achieve better vertical resolution in shallow water and areas with complex bathymetry.

### Dynamic Equations

The governing dynamical equations of three-dimensional, free-surface, Reynolds-averaged Navier-Stokes equations are (See appendix A for derivation and more details)

$$
\frac{\partial H_zu}{\partial t} + \frac{ \partial (uH_zu)}{\partial x} + \frac{ \partial( v H_z u)}{\partial y}  + \frac{\partial(\Omega H_zu)}{\partial \sigma}  - f H_z v = - \frac{ H_z}{\rho_0}\frac{ \partial p}{\partial x} -   H_z g\frac{ \partial \zeta }{\partial x} - \frac{\partial}{\partial \sigma}(-\frac{K_M}{H_z} \frac{\partial u}{\partial z} - \frac{\nu}{H_z}\frac{\partial u}{\partial \sigma})+F_u+D_u
$$

$$
\frac{\partial H_zv}{\partial t} + \frac{ \partial (uH_zv)}{\partial x} + \frac{ \partial( v H_z v)}{\partial y}  + \frac{\partial(\Omega H_zv)}{\partial \sigma}  + f H_z u = - \frac{ H_z}{\rho_0}\frac{ \partial p}{\partial y} -   H_z g\frac{ \partial \zeta }{\partial y} - \frac{\partial}{\partial \sigma}( -\frac{K_M}{H_z} \frac{\partial v}{\partial z} - \frac{\nu}{H_z}\frac{\partial v}{\partial \sigma})+F_v+D_v
$$

with the continuity equation

$$
\frac{\partial \zeta}{\partial t} + \frac{ \partial (H_zu)}{\partial x} + \frac{ \partial( H_z v)}{\partial y}  + \frac{\partial(H_z \Omega)}{\partial \sigma}  = 0
$$

equation of State

$$
\rho = \rho (T, S, P)
$$
and scalar transport equation for temperature and salinity

$$
\frac{\partial (H_z C)}{\partial t} + \frac{ \partial (u H_z C)}{\partial x} + \frac{ \partial( v H_z C)}{\partial y}  + \frac{\partial(\Omega H_z C)}{\partial \sigma} = -\frac{\partial}{\partial \sigma}( -\frac{K_C}{H_z} \frac{\partial C}{\partial z} - \frac{\nu_\theta}{H_z}\frac{\partial C}{\partial \sigma})
$$

Here, all variables are:


| name                                | Description                             |
| ----------------------------------- | --------------------------------------- |
| $u$                                 | horizontal velocity in x direction      |
| $v$                                 | horizontal velocity in y direction      |
| $\sigma$                            | the scaled sigma coordinate             |
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
#### Stretching Functions
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

### Turbulence Closure
When the Reyhnolds-averaged Navier-Stokes equations are first derived, there are terms including the average in the perturbations of velocity, such as $\overline{u'w'}$, making the number of unknown variables greater than the number of equations. In order to solve the problem, the turbulence closure technique is used, and the equations are:
$$
\overline{u'w'} = -K_M\frac{\partial u}{\partial z}\\
\overline{v'w'} = -K_M\frac{\partial v}{\partial z}\\
\overline{C'w'} = -K_C\frac{\partial C}{\partial z}
$$

As shown above, these equations approximate those terms with the gradient of velocity/tracer, implying that they flow down the local gradient of u, v, C respectively. This method of turbulence closure is also called the Gradient Transport Theory or K-theory.

## Data Assimilation with ROMS

### Why Data Assimilation is Needed?

Study of the data assimilation with partial observation is necessary and a fundamental challenge meteorology and oceanography because, in practice, it is impossible to measure the exact variable states of the entire system, due to but not limited to the following reasons:

1. Making measurements will cost too much. For example, measuring the velocity fields in the deep ocean.

2. The measurements have errors because of the limits of equipments.

3. The model may be inaccurate, such as making inappropriate assumptions.

There are several different data assimilation method could be applied in the ROMS. Here the two commonly used ones are simple nudging method and incremental 4D-VAR (I4D-VAR).

### Simple Nudging

The equation for simple nudging method is as follows.


$$
\frac{dx_a}{dt} = F_a(x(t))+g_l(t)(y(t)-x(t))\delta_{al}
$$

In the equations, subscripts a and l mean all variables and unobserved variables respectively. y is the observed data, and $g_l$ is the nudging coefficient. $\delta_{al}$ implies that we are only nudging the observed variables.

In this method, the nudging coefficient $g_l$ will modify the nudging strength and thus control the conditional Lyapunov exponent. If the largest Lyapunov exponent is smaller than zero in the system, all the unstable dimensions are constrained and, as a result, all the variables will be "nudged" to the right trojectory.

In order to constrain all the unstable dimensions, a minimum percentage of data is required to be observed. This is one of the question we would like to solve with ROMS: How many variables do we need to measure, to get a good prediction, or in other word, to constrain all the unstable dimensions?

However, when applying the nudging method, we may violate a specific physics law, since we are adding an extra nudging term g(y-x) to the system. For instance, if we are nudging the sea surface height (Choose $x_0$ to be $\zeta$), then by adding the extra term, we are violating the conservation of mass. Therefore, in order to satisfy the physics law, the dynamical nudging method, which we are planning to apply to ROMS, needs to be used. This method will be explained later.

### i-4DVar

TBD.

## Our Plan on Improving Data Assimilation in ROMS

### Time Delayed Nudging

In Zhe's paper [2], it was found that, by using the time delayed nudging method, the number of observed data required to make good predictions are reduced from 70% to 33% compared to the standard nudging method in the shallow water environment. It indicates that by introducing time delayed nudging into ROMS, we could potentially improve the system and make better forcast.

Adding explanation and equations for time delayed nudging (It'll be here in the next report...)

### Dynamical State and Parameter Estimation

Estimating parameters and unobserved state variables in nonlinear dynamical system is an essential aspect of the subject, and also a matter of interest to many other fields such as control theory, biological science and engineering [1]. In a common setting one has an experimental system described by a state vector $\vec{X}(t)$, which usually has a large dimensionality. However, it is common that only a sparse subset of $\vec{X}(t)$ could be recorded over time. For example, in the context of an ocean model, variables such as pressure and temperature can be easily probed on the surface of the ocean, and may be to a depth that is not too deep. This left us a big chanllenge in estimating the remaining dimensions of our state variable $\vec{X}(t)$. Once we have established a physical model for our system of interest, we also need to estimate, given the sparsely distributed experimental data, our model parameters $\vec{p}$.

Furthermore, if the model and the experimental system are chaotic, even if we have synchronized the data with our physical model, we would still face the problem that small perturbations in parameters or state variables can lead to large excursions near the synchronization manifold, and thus produce a very poor prediction to the future states.

An approach called dynamical state and parameter estimation addresses these instabilities and regularizes them, allowing for smooth surfaces in the space of parameters and initial conditions.

Let's say we have a model with dynamical variables $\vec{y}(t)$, and from experiments we have measured a subset of a dynamical variable $\vec{X}(t)$ from some initial time $t_I$ to a final time $t_F$. In order to make any prediction, one must estimate any unkown fixed parameters in the model as well as all the state variales at time $t_F$: $\vec{y}(t_F)$. Then the model-based predictions or forecasts for $t>t_F$ can be accomplished.

For simplicity we assume only the first dimension of $\vec{X}(t)$, i.e. $x_1(t)$ is measured over the time series. Also, we recognize that measurements are not made in continuous time but at discrete times. Thus we rewrite our dynamical state variable in discrete time: $\vec{X}(n)=\{x_1(n);\vec{X}_R(n)\}$, where $\vec{X}_R(n)$ are the unobserved dimensions.

Similarly, We wirte our D-dimensional model variable as $\vec{y}(n)=\{y_1(n);\vec{y}_R(n)\}$, where $y_1(n)$ corresponds to the observed $x_1(n)$, and the "rest" of them are indicated collectively by a subscript R. Last, we assume this model contains $L$ unknown parameters $\vec{p}=\{p_1,p_2,...,p_L\}$.

(Other details to be added here...)

The core equations for DSPE are:

$$
C(\vec{y};\vec{p};\vec{u})=\frac1{2N}\sum_{n=1}^{N-1} \{[x_1(n)-y_1(n)]^2+u^2(n)\}
y_1(n+1)=F_1(\vec{y}(n);\vec{p})+u(n)[x_1(n)-y_1(n)]
\vec{y}_R(n+1)=\vec{F}_R(\vec{y}(n);\vec{p})
$$

Where ${u}(n)$ is the regularization term. $F_1$ and $F_R$ are the dynamical equations of the system.

(To be completed):
1. How can we incorparate DSPE into ROMS.
2. What physics questions can be asked.


## Appendix

### A. Derivation of Reynolds-averaged Navier-Stokes equations

First, we shall start from the general Navier-Stokes equations and continuity equation in tensor notation,

$$
\frac{\partial}{\partial t}(\rho u_i)+u_j\frac{\partial}{\partial x_j}(\rho u_i)-\mu \frac{\partial^2}{\partial x_j\partial x_j}u_i = -\frac{\partial p}{\partial x_i}+f_i
$$
$$
\frac{\partial \rho}{\partial t}+\frac{\partial}{\partial x_i}(\rho u_i) = 0
$$

In the equations above, $\vec u$ is the velocity vector of the fluid and $u_i$ is the components in Cartesian coordinate. $\rho$, $\mu$, p are the density, dynamic viscosity, and pressure of the fluid respectively. $f_i$ is the external force term.

There are several approximations accepted by ROMS. First, the hydrostatic approximation is made, which implies that the horizontal scale is very large compared to the vertical scale. Furthermore, the Boussinesq approximation is assumed as well, which means that the the variation in density is only important in the buoyancy term. As a result, the equations can be simplified as

$$
\frac{\partial u_i}{\partial t}+u_j\frac{\partial u_i}{\partial x_j}-\nu \frac{\partial^2}{\partial x_j\partial x_j}u_i = -\frac{\rho}{\rho_0}g_i+f_i
$$
$$
\frac{\partial u_i}{\partial x_i} = 0
$$
$\nu = \mu/\rho_0$ is the kinematic viscosity.

Then, we need to apply Reynolds decomposition to the systems, namely that we will decompose $\vec u(x,y,z,t)$ and other variables as

$$
\vec u(x,y,z,t) = \overline{u(x,y,z,t)}+u'(x,y,z,t)
$$
$\overline{u(x,y,z,t)}$ is the average of the velocity, while $u'(x,y,z,t)$ is the fluctuating term, satisfying $\overline{u'(x,y,z,t)} = 0$. Then the equations will become

$$
\frac{\partial (\overline{u_i}+u_i')}{\partial t}+(\overline{u_j}+u_j')\frac{\partial (\overline{u_i}+u_i')}{\partial x_j}-\nu \frac{\partial^2}{\partial x_j\partial x_j}(\overline{u_i}+u_i') = -\frac{\rho}{\rho_0}g_i+(\overline{f_i}+f_i')
$$
$$
\frac{\partial (\overline{u_i}+u_i')}{\partial x_i} = 0
$$

Next, we will take the average of the equations, taking advantage of $\overline{u'(x,y,z,t)} = 0$, and rewriting $\overline{u_i}$ and $\overline{f_i}$ as $u_i$ and $f_i$ for simplicity, we have

$$
\frac{\partial u_i}{\partial t}+u_j\frac{\partial u_i}{\partial x_j}+\frac{\partial}{\partial x_j}\overline{u_i'u_j'}-\nu \frac{\partial^2}{\partial x_j\partial x_j}u_i = -\frac{\rho}{\rho_0}g_i+f_i
$$
$$
\frac{\partial u_i}{\partial x_i} = 0
$$

After applying hydrostatic approximation and turbulence closure (as discussed in the main text), the equations are
$$
\frac{\partial u}{\partial t}+\vec v\cdot\nabla u-fv =-\frac{1}{\rho_0}\frac{\partial p}{\partial x}+\frac{\partial^2 u}{\partial z^2}(\nu+K_M)+F_u+D_u
$$
$$
\frac{\partial u}{\partial t}+\vec v\cdot\nabla v+fu =-\frac{1}{\rho_0}\frac{\partial p}{\partial y}+\frac{\partial^2 v}{\partial z^2}(\nu+K_M)+F_v+D_v
$$
$$
\frac{\partial p}{\partial z} = -\rho g
$$
$$
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z} = 0
$$
In the equations above, $\vec v$ is the velocity vector $(u,v,w)$, while u, v, w are the x,y,z velocity components. -fv and +fu are the Coriolis force terms. F is the external forcing, and D is the optional extra diffusing term.

Last, we shall apply the terrain-following coordinate and introduce the stretching factor $H_z = \frac{\partial z}{\partial \sigma}$, as desctribed in "Vertical S-coordinate" section. Then using the chain rules,
$$
{(\frac{\partial}{\partial x})}_{z} = {(\frac{\partial}{\partial x})}_{\sigma}-(\frac{1}{H_z}){(\frac{\partial z}{\partial x})}_{\sigma}\frac{\partial}{\partial \sigma}
$$
$$
{(\frac{\partial}{\partial y})}_{z} = {(\frac{\partial}{\partial y})}_{\sigma}-(\frac{1}{H_z}){(\frac{\partial z}{\partial y})}_{\sigma}\frac{\partial}{\partial \sigma}
$$
$$
\frac{\partial}{\partial z} = \frac{1}{H_z}\frac{\partial}{\partial \sigma}
$$
we can get the equations in the "Dynamic Equations" section.

### B. Details of the I4D-VAR Code (Reference or starting point for us to apply our dynamical nudging or DSPE method to the ROMS)

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

Reference

[1] Abarbanel, Henry DI, et al. "Dynamical state and parameter estimation." SIAM Journal on Applied Dynamical Systems (2009): 1341-1381.
[2] An, Z., Rey, D., Ye, J., and Abarbanel, H. D. I.: Estimating the state of a geophysical system with sparse observations: time delay methods to achieve accurate initial states for prediction, Nonlin. Processes Geophys., 24, 9-22, doi:10.5194/npg-24-9-2017, 2017.
