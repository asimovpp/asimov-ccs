%% plot-structured.m
%
% Plots the PETSc binary output for a structured grid
% Expects U, V, P
%

%% Setup paths, load PETSc functions
petsc_path = "/work/Projects/ASiMoV/petsc"; % User-specific
petsc_mlab = [petsc_path "/share/petsc/matlab"];
%%addpath(petsc_mlab);
%%source([petsc_mlab "/PetscBinaryRead.m"]);

%% Set mesh size
cps = 129; % Cells per side (assumed square mesh)
L   = 1.0; % Length of side (assumed square mesh)
dx  = L / cps;
x   = linspace(0.5 * dx, L - 0.5 * dx, cps);
mp  = floor((cps + 1) / 2);

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

quiver(U(1:3:cps, 1:3:cps), V(1:3:cps, 1:3:cps), 5);
axis tight;
print(fig, "quiver.png", "-dpng")

velmag = sqrt(U.**2 + V.**2);
pcolor(velmag);
colorbar();
colormap("jet");
shading interp;
print("velmag", "-dpng");
contour(velmag);
print contour.png

%% Ghia data (Re=100)

% u - plotted through vertical centreline
yu = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1];
ug = [0 -0.03717 -0.04192 -0.04775 -0.06434 -0.10150 -0.15662 -0.21090 -0.20581 -0.13641 0.00332 0.23151 0.68717 0.73722 0.78871 0.84123 1];

xv = [0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5 0.8047 0.8594 0.9063 0.9453 0.9531 0.96099 0.9688 1];
vg = [0 0.09233 0.10091 0.10890 0.12317 0.16077 0.17507 0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0];

%% Ghia data (Re=1000)

% u - plotted through vertical centreline
% yu = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1];
% ug = [0 -0.18109 -0.20196 -0.22220 -0.29730 -0.38289 -0.27805 -0.10648 -0.06080 0.05702 0.18719 0.33304 0.46604 0.51117 0.57492 0.65928 1];

% xv = [0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5 0.8047 0.8594 0.9063 0.9453 0.9531 0.96099 0.9688 1];
% vg = [0 0.27485 0.29012 0.30353 0.32627 0.37095 0.33075 0.32235 0.02526 -0.31966 -0.42665 -0.51550 -0.39188 -0.33714 -0.27669 -0.21388 0];

%% Ghia data (Re=10000)

% u - plotted through vertical centreline
% yu = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1];
% ug = [0 -0.42735 -0.42537 -0.41657 -0.38000 -0.32709 -0.23186 -0.07540 0.03111 0.08344 0.20673 0.34635 0.47804 0.48070 0.47783 0.47221 1];

% xv = [0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5 0.8047 0.8594 0.9063 0.9453 0.9531 0.96099 0.9688 1];
% vg = [0 0.43983 0.43733 0.43124 0.41487 0.35070 0.28003 0.27224 0.00831 -0.30719 -0.36737 -0.41496 -0.45863 -0.49099 -0.52987 -0.54302 0];

%% Plot against Ghia
plot(x, U(:, mp));
hold on;
plot(yu, ug);
print("u-line", "-dpng");
hold off

plot(x, V(mp, :));
hold on;
plot(xv, vg);
print("v-line", "-dpng");
