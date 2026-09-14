---
title: Transient discretisation schemes
---


Transient schemes are implemented using the kernel approach. 
The `transient_kernels` module defines an abstract type `transient_kernel` for which the initialisation function `init` is defined for each specific scheme. This initialisation sets private quantities defining the different coefficients (implicit, explicit) and variables (theoretical order and width) that are unique to each scheme.

Kernels that require more than a single previous time step to populate their stencil need to revert back to a lower order scheme (requiring a smaller stencil) for the very begining of the simulation. This is done as part of the initialisation setup by having the implicit and explicit coefficients initialised by 1D and 2D arrays respectively. For example `implicit_coeffs_trans(i)` contains the implicit coefficient for the scheme used for the i-th timestep, with the last element of this array being the coefficient of the scheme itself.



# Discretisation schemes

## 1st order (Implicit Euler)

Implicit Euler is a first order scheme that uses a single older timestep. Following the series expansion:

$$
\frac{\phi^{n+1} - \phi^{n}}{\Delta t} = \frac{\partial \phi}{\partial t} + O(\Delta t)
$$

This leads to a trivial approximation of the derivative. Implicit and explicit coefficients are then defined as:
```Fortran
implicit_coeff = 1.0
explicit_coeffs(1) = 1.0
```

## 2nd order (Three time level)

Implicit three time level is a second order scheme requiring two older timesteps. It uses a quadratic backward approximation in time:

$$
\frac{\phi^{n+1} - \phi^{n}}{\Delta t} = \frac{\partial \phi}{\partial t} -  \frac{\partial^2 \phi}{\partial t^2}\frac{\Delta t}{2} + O(\Delta t^2)
$$

Considering the same approximation between time steps $n+1$ and $n-1$ (separated by $2\Delta t$) and substracting the expressions leads to

$$
\frac{\partial^2 \phi}{\partial t^2} = \frac{\phi^{n+1} - 2\phi^n + \phi^{n-1}}{\Delta t^2} + O(\Delta t^2)
$$
Which can be substituted in the first expression leading to:

$$
\frac{\partial \phi}{\partial t} = \frac{3\phi^{n+1} -4\phi^n + \phi^{n-1}}{2\Delta t} + O(\Delta t^2)
$$

This leads to an implicit coefficient of $\frac{3}{2}$ and explicit coefficients of $2$ and $-\frac{1}{2}$, defined in ccs as:
```Fortran
implicit_coeff = 1.5
explicit_coeffs(1) = 2.0
explicit_coeffs(2) = -0.5
```

## Theta scheme

This scheme is a blend between the first order Implicit Euler scheme and the second order three time level method. The blending factor, $\theta$, can vary between 0 and 1, 0 leading to the implicit Euler scheme and 1 leading to the three time level scheme:

$$
\frac{\partial \phi}{\partial t} \rightarrow \frac{(2 + \theta)\phi^{n+1} -2(1+\theta)\phi^n + \theta\phi^{n-1}}{2\Delta t} 
$$
The order of the scheme when $\theta$ is strictly between 0 and 1 isn't clearly defined here but is at least 1. In ccs implicit and explicit coefficients are defined as follow:

```Fortran
implicit_coeff = 1 + 0.5*theta
explicit_coeffs(1) = 1.0 + theta
explicit_coeffs(2) = -0.5*theta
```