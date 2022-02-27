function [u]=impl_BEC_linearized(nx, ny, nz, xl, xr, yl, yr, zl, zr, eigvec, lam, t_lambda)

[nodes, elements] = initialize(nx, ny, nz, xl, xr, yl, yr, zl, zr);

rho = eigvec.^2;
[mat_AA_lap,mat_sec_AA,~] =eig_gstif1(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements);
mat_AA = eig_gstif2(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements,rho);

% f = lam*mat_sec_AA*eigvec - mat_AA*eigvec;
% u = cg_solve(mat_AA_lap,f);
f = (lam+t_lambda)*mat_sec_AA*eigvec - mat_AA*eigvec;
u = cg_solve(mat_AA_lap+t_lambda*mat_sec_AA, f);

return