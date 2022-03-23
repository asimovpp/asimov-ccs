%% plot-structured.m
%
% Plots the PETSc binary output for a structured grid
% Expects U, V, P
%

%% Setup paths, load PETSc functions
petsc_path = "~/.local/opt/petsc/v3.16.5"; % User-specific
petsc_mlab = [petsc_path "/share/petsc/matlab"];
addpath(petsc_mlab);
source([petsc_mlab "/PetscBinaryRead.m"]);

%% Set mesh size
cps = 129; % Cells per side (assumed square mesh)
L   = 1.0; % Length of side (assumed square mesh)

%% Load data
p = PetscBinaryRead("p");
u = PetscBinaryRead("u");
v = PetscBinaryRead("v");

%% Reshape data
P = reshape(p, cps, cps)';
U = reshape(u, cps, cps)';
V = reshape(v, cps, cps)';

%% Plot
fig = figure();
set(fig, "Visible", "Off");

pcolor(P);
colorbar();
colormap("jet");
shading interp;
print(fig, "p.png", "-dpng")

pcolor(U);
colorbar();
colormap("jet");
shading interp;
print(fig, "u.png", "-dpng")

pcolor(V);
colorbar();
colormap("jet");
shading interp;
print(fig, "v.png", "-dpng")

quiver(U(1:2:cps, 1:2:cps), V(1:2:cps, 1:2:cps), 5);
axis tight;
print(fig, "quiver.png", "-dpng")
