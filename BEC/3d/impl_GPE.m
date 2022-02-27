function [mat_stif,mat_sec_AA,mat_H,eigval0,eigvec] = impl_GPE(nx,ny,nz,xl,xr,yl,yr,zl,zr,t_lambda)
% To solve GPE for BEC by standard FEM
% Fang Liu, 2021/02/22

% clear;
% t_lambda=50.0;
% 
% nx = 10;
% ny = nx;
% nz = nx;
% l = -5.0;
% r = -l;
% 
% xl =l;% 
% xr =r;% 
% yl =l;% 
% yr =r;%
% zl =l;% 
% zr =r;% the bound of the computational domain

[nodes, elements] = initialize(nx, ny, nz, xl, xr, yl, yr, zl, zr);

dof=(nx-1)*(ny-1)*(nz-1);
nne=(nx+1)*(ny+1)*(nz+1);

err_eigenvalue=1;
rho = zeros(dof,1);
rho0 = zeros(dof,1);

Fm = zeros(dof,1);
Fm0 = zeros(dof,1);

eigval0=0.0;
self_iter_time=0;
self_iter_max=200;
err_max=1.e-5;

beta=200.0;
shell=1.0;

uu=zeros(self_iter_max,dof);
vv=zeros(self_iter_max,dof);

while (self_iter_time<self_iter_max) && (err_eigenvalue>err_max)
 
    if self_iter_time==0
        [mat_AA_lap,mat_sec_AA,mat_H] =eig_gstif1(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements);
    end
    mat_AA = eig_gstif2(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements,rho);
            
    mat_stif= mat_AA_lap + mat_AA + t_lambda*mat_sec_AA;
    %condest(mat_stif)
    
    [eigvec,eigval] = eigs(mat_stif,mat_sec_AA,1,'SM'); 
    eigval=eigval-t_lambda;       
    err_eigenvalue=abs(eigval-eigval0);
    eigval0=eigval;
    
    fprintf('%i-th iteration, err_eigval = %f, eigenvalue = %f\n', self_iter_time, err_eigenvalue, eigval);
    
    mod_u = module_u(nx, ny, nz, xl, xr, yl, yr, zl, zr, eigvec);
    eigvec=eigvec/mod_u;
    
% method 1     
    if (err_eigenvalue > err_max && self_iter_time <(self_iter_max-1))
        [Fm,vv,uu,rho,rho0]=new_rho(err_eigenvalue, dof, Fm, eigvec,self_iter_time,self_iter_max,rho,rho0,uu,vv,shell,err_max);
    end
    
    self_iter_time=self_iter_time+1;
    
end

E = energy(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements, eigvec);


    
    
