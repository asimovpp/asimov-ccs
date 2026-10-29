---
title: Residual monitoring and storage
---

Residual monitoring and field output serve different purposes. Monitoring reduces
cell residuals to scalar convergence measures. Setting `variables.store_residuals`
to `true` additionally retains residual fields for output; it must not change the
equations solved, their right-hand sides, or the computed solution.

For a matrix equation, the stored residual is the raw cell vector

\[
  r = b - A\phi.
\]

The SIMPLE momentum residual is evaluated before the linear solve, so that it
measures the current nonlinear iterate. Transported scalar residuals are evaluated
after their linear solve. Scalar-only runs can store these residuals without a
nonlinear convergence monitor.

Pressure correction solves an error equation. The cell mass imbalance is therefore
the source of its right-hand side. The pressure residual field records this raw
imbalance before pressure-correction preparation removes the mean, where required,
and applies the right-hand-side sign convention. It is an independent snapshot:
enabling output does not redirect assembly away from the solver's right-hand side.

After correcting face fluxes, SIMPLE also monitors the corrected cell mass
imbalance under `mf`. This is a cell-based diagnostic, whereas `mf` itself is
face-based. Only its norms are recorded; the cell vector is not copied into the
face residual storage.

For monitoring, each owned-cell residual is divided by the cell volume to the
power of two thirds. The quantity labelled `L2` in the residual log is the global
root mean square of those normalised values. `Linf` is their maximum absolute
value. Halo cells are excluded from both reductions. The retained output fields
contain raw residuals rather than these volume-normalised values.

The residuals module owns temporary vector storage for monitoring when no output
field is allocated. It creates that workspace lazily with the equation vector's
parallel layout and reuses it within the nonlinear solve. Initialisation and
finalisation release the workspace, so it cannot retain a previous mesh or
communicator across solver lifetimes. A storage-only scalar calculation writes
directly into the field and needs no temporary residual vector.
