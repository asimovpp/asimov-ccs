---
title: Source term treatment
---

In ASiMoV-CCS we are concerned with solving transport equations (momentum, scalar transport) that
are of the general form
\begin{equation}
	\frac{\partial\phi}{\partial{}t} + \nabla\cdot\phi{}u = \nabla\cdot\Gamma\nabla\phi + S
\end{equation}
where \(\phi\) is the transported variable, u the velocity vector, \Gamma the diffusion coefficient
and S is a source term.
Frequently source terms have linear and fixed components: \(S=S^{\phi}-R\phi{}\) so that the equation
is written as
\begin{equation}
	\frac{\partial\phi}{\partial{}t} + \nabla\cdot\phi{}u + R\phi = \nabla\cdot\Gamma\nabla\phi + S^{\phi}
\end{equation}
which can be beneficial to the discretisation.

Following a standard finite volume method we obtain a discretised equation of the form
\begin{equation}
	\frac{\partial\phi_P}{\partial{}t}V_P + \sum_{f}\phi{}u_{f}\cdot{}n_{f}A_f + R_{P}\phi_{P}V_P = 
		\sum_{f}\Gamma\nabla\phi_{f}n_{f}A_f + S^{\phi}_{P}V_P
\end{equation}
which yields the linear system
\begin{equation}
	(A + RV)\phi = (b + SV)
\end{equation}
where the formation of \(A\) and \(b\) have been discussed elsewhere.
The additional terms \(RV\) and \(SV\) are a diagonal matrix and vector, respectively, representing
the *integrated* source terms; this reflects the code implementation which expects integrated
sources, *i.e.* user and model code evaluating source terms should return \(RV\) and \(SV\).
This is demonstrated in the body of `zero_sources`:
```
do index_p = 1, local_num_cells
  call create_cell_locator(index_p, loc_p)
  call get_volume(loc_p, V_p)
  R_data(index_p) = 0 * V_p
  S_data(index_p) = 0 * V_p
end do
```
even though this subroutine is provided to simplify the case when no source is required, it
multiplies the zero point-wise value by the cell volume to serve as a template for other cases.
