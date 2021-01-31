% u = [0.0022, 5.8540e-04,  1.4866e-04, 3.7344e-05 ];

% h= 1./[8, 16, 32, 64];
% u = [ 0.0094, 0.0030, 8.2173e-04, 2.1225e-4];
% u = [0.0017, 4.36e-4, 1.09e-4];
% u = [0.0026, 8.4213e-04, 2.3115e-04, 5.9254e-5];
% ep = [0.0043, 0.0012, 3.2833e-4, 8.333e-5];
% figure()
% loglog(h, ep, '-o')
% hold on
% loglog(h, 0.2*h.^1.9, 'k')
% shg

%% 
% Pde=p;
load ('mesh_store_m8.mat')
load pde_store_m8.mat
mesh8=Meshes;
pde8=Pde;
load ('mesh_store_m16.mat')
load pde_store_m16.mat
mesh16=Meshes;
pde16=Pde;
load ('mesh_store_m32.mat')
load pde_store_m32.mat
mesh32=Meshes;
pde32=Pde;
load ('mesh_store_m64.mat')
load pde_store_m64.mat
mesh64=Meshes;
pde64=Pde;
load pde_store_m128.mat;
load mesh_store_m128.mat;
pde128=Pde;
mesh128=Meshes;

%% Interpolation

pde8.U_interp = interp2(full(mesh128.X1_U), full(mesh128.Y1_U), full(pde128.U1), full(mesh8.X1_U), full(mesh8.Y1_U), 'cubic');
pde16.U_interp = interp2(full(mesh128.X1_U), full(mesh128.Y1_U), full(pde128.U1), full(mesh16.X1_U), full(mesh16.Y1_U), 'cubic');
pde32.U_interp = interp2(full(mesh128.X1_U), full(mesh128.Y1_U), full(pde128.U1), full(mesh32.X1_U), full(mesh32.Y1_U), 'cubic');
pde64.U_interp = interp2(full(mesh128.X1_U), full(mesh128.Y1_U), full(pde128.U1), full(mesh64.X1_U), full(mesh64.Y1_U), 'cubic');

%% Computing Error
dx=1/8;
dy=1/8;
err_u8 = sqrt(sum(sum((pde8.U1-pde8.U_interp).^2))*dx*dy)

dx=1/16;
dy=1/16;
err_u16 = sqrt(sum(sum((pde16.U1-pde16.U_interp).^2))*dx*dy)

dx=1/32;
dy=1/32;
err_u32 = sqrt(sum(sum((pde32.U1-pde32.U_interp).^2))*dx*dy)

dx=1/64;
dy=1/64;
err_u64 = sqrt(sum(sum((pde64.U1-pde64.U_interp).^2))*dx*dy)

%% Plot

h = 1./[8, 16, 32, 64];
u = [err_u8, err_u16, err_u32, err_u64];
figure
loglog(h, u, '-o')
hold on
loglog(h, h.^2, 'k')
hold off
shg








