function [u]=impl_GPE_linearized(nx, ny, xl, xr, yl, yr, eigvec, lam, t_lambda)

[nodes, elements] = initialize(nx, ny, xl, xr, yl, yr);

rho = eigvec.^2;
[mat_AA_lap,mat_sec_AA,~] =eig_gstif0(nx, ny, xl, xr, yl, yr, nodes, elements);
mat_AA = eig_gstif2(nx, ny, xl, xr, yl, yr, nodes, elements,rho);

% f = lam*mat_sec_AA*v - mat_AA*v;
% u = cg_solve(mat_AA_lap,f);

f = (lam+t_lambda)*mat_sec_AA*eigvec - mat_AA*eigvec;
u = cg_solve(mat_AA_lap+t_lambda*mat_sec_AA, f);

return