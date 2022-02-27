function [mat_stif,mat_sec_AA,mat_H,eigval0,eigvec] = impl_oscillator(nx, ny, nz, xl, xr, yl, yr, zl, zr,t_lambda)
% To solve equation for oscillator by standard FEM
% Fang Liu, 2021/02/19

% clear;
% t_lambda=4.0;
% 
% nx = 20;
% ny = nx;
% nz = nx;
% l = -5;
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

eigval0=0.0;
self_iter_time=0;
self_iter_max=100;
err_max=1.e-6;

while (self_iter_time<self_iter_max) && (err_eigenvalue>err_max)
 
    if self_iter_time==0
        [mat_AA_lap,mat_sec_AA,mat_H] =eig_gstif0(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements);
    end
    mat_AA = eig_gstif2(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements,rho);
            
    mat_stif = mat_AA_lap + mat_AA + t_lambda*mat_sec_AA;
    
    [eigvec,eigval] = eigs(mat_stif,mat_sec_AA,1,'SM');
%     [eigvec,eigval] = inv_power(mat_stif,mat_sec_AA);
    eigval=eigval-t_lambda;
    
    err_eigenvalue=abs(eigval-eigval0);
    eigval0=eigval;
    
    fprintf('%i-th iteration, err_eigval = %f, eigenvalue = %f\n', self_iter_time, err_eigenvalue, eigval);
    
    mod_u = module_u(nx, ny, nz, xl, xr, yl, yr, zl, zr, eigvec);
    eigvec=eigvec/mod_u;
    
% method 1     
%     if (err_eigenvalue >err_max && self_iter_time <(self_iter_max-1))
%         [Fm,vv,uu,rho,rho0]=new_rho(err_eigenvalue, dof, Fm, eigvec,self_iter_time,self_iter_max,rho,rho0,uu,vv);
%     end

% % method 2
    q = 1;
    rho_out = eigvec.^2;
    rho = q*rho_out + (1-q)*rho;
    
    self_iter_time=self_iter_time+1;
    
end



    
    
