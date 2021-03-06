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

### Boundary conditions

#### Vertical boundary conditions

In ROMS, there is a bottom layer is assume to have zero velocity (no-slip boundary condition), the thickness of the layer is defined as the roughness of the bottom $z_0$. However, usually $z_0$ is much smaller than the thickness of one vertical layer in the discretized vertical s-coordinate and the horizontal velocity (u,v) are evaluated at the mid point in each vertical layer ($\rho$ point in the following graph). For example, the $z_0$ value used in the previous toy model is 20 cm, while the thickness of the bottom layer is around 1000m.  Then we need to estimate the effect of the bottom no-slip layer on the lowest $\rho$ points, which gives us the bottom boundary condition to use in u and v.

![Vertical grid](1.png)

The method in ROMS to estimate the bottom stress from the no-slip boundary condition is to assume a layer of constant Reynolds stress near the bottom [3]. Then applying the turbulence closure we discussed before, we have

$$
K_M\frac{\partial u}{\partial z} = \tau_b^x(x,y,t)
$$
$$
K_M\frac{\partial v}{\partial z} = \tau_b^y(x,y,t)
$$

Next, in order to get the bottom stress $\tau_b$, we assumed a logarithmic velocity profile of the bottom velocity, which is an analogy to the wind stress effect at the surface layer. The profile satisfies

$$
v(z) = \frac{v_\star}{\kappa}log(\frac{z}{z_0})
$$

In the equation above, z is the height above the bottom, which is ($z_{\rho 1}-z_{w0}$) for our case. $z_0$ is the roughness of the bottom as mentioned above. $\kappa$ is the von Karman's constant, and $v_\star$ is the current friction velocity defined as $\rho v_\star^2 = \tau_b$. Therefore, since we have the initial velocity profile, we can use the logarithmic profile above to obtain the estimate of bottom stress, then apply it to the bottom boundary conditions, which is actually what ROMS does in the code.

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

###Twin experiments

The twin experiments we set up is the same as our previous work, except that the wind force is added. The initial conditions are shown below. one set start with an initial temperature gradient and 0 initial velocity, while the other one has a Gaussian noise added on temperature and velocity.

<img src="Ti_Data.png" alt="Drawing" style="width: 400px;"/><img src="Ti_Model.png" alt="Drawing" style="width: 400px;"/>
<img src="U_Data.png" alt="Drawing" style="width: 400px;"/><img src="U_Model.png" alt="Drawing" style="width: 400px;"/>
<img src="V_Data.png" alt="Drawing" style="width: 400px;"/><img src="V_Model.png" alt="Drawing" style="width: 400px;"/>

The wind force added has the form $\tau_u (i,j) = -0.1cos(2\pi y(i,j)/L)\space m^2/s^2$, where i,j are the index of grid points in x and y directions respectively, L is the total length in y direction, and subscript u indicates that the wind force is in x direction. The boundary conditions are periodical for y direction and closed for x direction. The bottom is flat. The other relevant parameters are listed below, as well as a plot for the wind stress.

| name                                | Description                             |Value             |
| ----------------------------------- | --------------------------------------- |------------------|
| $N_i$                               | Number of x direction $\rho$ points     |200               |
| $N_j$                               | Number of y direction $\rho$ points     |100               |
| $N_\sigma$                          | Number of vertical layers               |20                |
| $dt$                                | Time step size                          |600s (10 min)     |
| $N_{time}$                          | Number of time steps                    |42200             |
| $N_{his}$                           | Number of time steps between observation|144 (1 day)       |
| $Zo_b$                              | Bottom Roughness                        |0.02m             |
| $\theta_s$                          | See Vertical S-coordinate section       |7                 |
| $\theta_b$                          | See Vertical S-coordinate section       |0.1               |

<img src="Windstress.png" alt="Drawing" style="width: 600px;"/>

\
The total number of unknown variables are: 200x100x20x3(u, v, and temperature)+200x100($\zeta$) ~
 24.4 million. Below the graphs for the final velocity and temperature are provided.

 <img src="Tf_Data.png" alt="Drawing" style="width: 400px;"/><img src="Tf_Model.png" alt="Drawing" style="width: 400px;"/>
 <img src="Uf_Data.png" alt="Drawing" style="width: 400px;"/><img src="Uf_Model.png" alt="Drawing" style="width: 400px;"/>
 <img src="Vf_Data.png" alt="Drawing" style="width: 400px;"/><img src="Vf_Model.png" alt="Drawing" style="width: 400px;"/>


###Simple Nudging
(We will be doing this for the next week, and here's the plan)

We would like to apply the simple nudging method to the twin experiments, and see how many percentage of data we need to get good predictions after nudging. Right now based on the previous results in our group for shallow water equations, I would expect something around 70% to 80% percent.
However, We are thinking on how to pick up the points, since the data in deep ocean is usually hard to get. We are considering only do nudging on the layers from top to down. For example, if we use 80% of the data, then we can use the top 18 layers, which is 90%, and then pick 90% of the points randomly from the horizontal plane, then we will have 90%x90%=81% observed data, which is close to 80.)

### Time Delayed Nudging

In Zhe's paper [2], it was found that, by using the time delayed nudging method, the number of observed data required to make good predictions are reduced from 70% to 33% compared to the standard nudging method in the shallow water environment. It indicates that by introducing time delayed nudging into ROMS, we could potentially improve the system and make better forcast.

Adding explanation and equations for time delayed nudging...

### Dynamical State and Parameter Estimation

Estimating parameters and unobserved state variables in nonlinear dynamical system is an essential aspect of the subject, and also a matter of interest to many other fields such as control theory, biological science and engineering [1]. In a common setting one has an experimental system described by a state vector $\vec{X}(t)$, which usually has a large dimensionality. However, it is common that only a sparse subset of $\vec{X}(t)$ could be recorded over time. For example, in the context of an ocean model, variables such as pressure and temperature can be easily probed on the surface of the ocean, and may be to a depth that is not too deep. This left us a big chanllenge in estimating the remaining dimensions of our state variable $\vec{X}(t)$. Once we have established a physical model for our system of interest, we also need to estimate, given the sparsely distributed experimental data, our model parameters $\vec{p}$.

Furthermore, if the model and the experimental system are chaotic, even if we have synchronized the data with our physical model, we would still face the problem that small perturbations in parameters or state variables can lead to large excursions near the synchronization manifold, and thus produce a very poor prediction to the future states.

An approach called dynamical state and parameter estimation addresses these instabilities and regularizes them, allowing for smooth surfaces in the space of parameters and initial conditions.

Let's say we have a model with dynamical variables $\vec{y}(t)$, and from experiments we have measured a subset of a dynamical variable $\vec{X}(t)$ from some initial time $t_I$ to a final time $t_F$. In order to make any prediction, one must estimate any unkown fixed parameters in the model as well as all the state variales at time $t_F$: $\vec{y}(t_F)$. Then the model-based predictions or forecasts for $t>t_F$ can be accomplished.

For simplicity we assume only the first dimension of $\vec{X}(t)$, i.e. $x_1(t)$ is measured over the time series. Also, we recognize that measurements are not made in continuous time but at discrete times. Thus we rewrite our dynamical state variable in discrete time: $\vec{X}(n)=\{x_1(n);\vec{X}_R(n)\}$, where $\vec{X}_R(n)$ are the unobserved dimensions.

Similarly, We write our D-dimensional model variable as $\vec{y}(n)=\{y_1(n);\vec{y}_R(n)\}$, where $y_1(n)$ corresponds to the observed $x_1(n)$, and the "rest" of them are indicated collectively by a subscript R. Last, we assume this model contains $L$ unknown parameters $\vec{p}=\{p_1,p_2,...,p_L\}$.

In discrete time limit, we can represent our physical model as a set of equations as a function of discrete time:

$$
y_1(n+1)=F_1(\vec{y}(n);\vec{p})
$$
$$
\vec{y}_R(n+1)=\vec{F}_R[\vec{y}(n);\vec{p}]
$$

We see from above that the phase-space behavior of our model is determined by initial state $\vec{y}(0)$ and model parameters $\vec{p}=\{p_1,p_2,...,p_L\}$.

Recall that our goal is to synchronize the (rather limited) measurements with our model, to the best of our knowledge, a standard way to quantify this is to introduce the Mean-Squared Error Function:

$$
C(\vec{y}(0);\vec{p})\equiv\frac1{2N}\sum_{n=1}^{N-1} \{[x_1(n)-y_1(n)]^2\}
$$

This is a function of all the parameters $\vec{p}$ because the time evolution of our model depends on these parameters, and this is a function of $\vec{y}_R(0)$ because the time evolution requires full initial states. With these said, it seem that if we can tune our model parameter and initial state to an optimal setting that minimizes the Error Function $C$, we will be done the estimation process. However, in real world, it is rare that a physical system can be this perfect -- i.e., usually these systems are chaotic, generating wild deviations between measurements and predictions. Mathematically speaking, it is because the synchronization manifold $\vec{X}(n)\simeq\vec{y}(n)$ is not stable. This indicates, on the other hand, at least one positive Lyapunov exponent in Lyapunov spectrum.

This drive us to investigate the stability and regularization of the synchronization manifold, which is a must have for minimization of $C(\vec{y}(0);\vec{p})$

(Mar 3rd: further details to be added here... Stabilization, How DSPE reduces CLE, details of the DSPE equations etc.)

The core equations for DSPE are:

$$
C(\vec{y};\vec{p};\vec{u})=\frac1{2N}\sum_{n=1}^{N-1} \{[x_1(n)-y_1(n)]^2+u^2(n)\}
$$
$$
y_1(n+1)=F_1(\vec{y}(n);\vec{p})+u(n)[x_1(n)-y_1(n)]
$$
$$
\vec{y}_R(n+1)=\vec{F}_R[\vec{y}(n);\vec{p}]
$$

Where $C$ is the cost function we need to minimize, $u(n)$ is the regularization term. $F_1$ and $F_R$ are the dynamical equations of the system.

(To be completed Mar 3rd):
1. How can we incorparate DSPE into ROMS.
2. What physics questions can be asked.
3. DSPE with time-delayed nudging.

### Discussion: the Number of Required Observations in Data Assimilation in ROMS

An important Chapter to be started in Mar. 3rd report. Reference: [4]

[Simple nudging vs. DSPE+time-delayed nudging]

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



Reference

[1] Abarbanel, Henry DI, et al. "Dynamical state and parameter estimation." SIAM Journal on Applied Dynamical Systems (2009): 1341-1381.
[2] An, Z., Rey, D., Ye, J., and Abarbanel, H. D. I.: Estimating the state of a geophysical system with sparse observations: time delay methods to achieve accurate initial states for prediction, Nonlin. Processes Geophys., 24, 9-22, doi:10.5194/npg-24-9-2017, 2017.
[3] Dewey, Richard K., and William R. Crawford. "Bottom stress estimates from vertical dissipation rate profiles on the continental shelf." Journal of Physical Oceanography 18.8 (1988): 1167-1177.
[4] Whartenby, William G., John C. Quinn, and Henry DI Abarbanel. "The number of required observations in data assimilation for a shallow-water flow." Monthly Weather Review 141.7 (2013): 2502-2518.
