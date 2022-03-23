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
pcolor(P);
colorbar();
colormap("jet");
shading interp;
print p.png

pcolor(U);
colorbar();
colormap("jet");
shading interp;
print u.png;

pcolor(V);
colorbar();
colormap("jet");
shading interp;
print v.png;

quiver(U(1:2:cps, 1:2:cps), V(1:2:cps, 1:2:cps), 5);
axis tight;
print quiver.png
