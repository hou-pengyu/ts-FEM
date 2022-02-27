function [u]=impl_he_atom_linearized(nx, ny, nz, xl, xr, yl, yr, zl, zr, eigvec, eigval, t_lambda)

[nodes, elements] = initialize(nx, ny, nz, xl, xr, yl, yr, zl, zr);

rho = 2*eigvec.^2;
[mat_AA_lap,mat_AA_lap1,mat_sec_AA] = eig_gstif0(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements);
V_Har = V_Har_pot(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements, rho, mat_AA_lap);
mat_AA = eig_gstif2(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements, V_Har, rho);

% f = (eigval)*mat_sec_AA*eigvec - mat_AA*eigvec;
% u = cg_solve(mat_AA_lap1,f);
%Fang
f = (eigval+t_lambda)*mat_sec_AA*eigvec - mat_AA*eigvec;
u = cg_solve((mat_AA_lap1+t_lambda*mat_sec_AA),f);

return