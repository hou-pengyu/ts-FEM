function [f_m,v,u,rho_m1]=broyden(m, f_m,f_m0, v, u, rho_m,rho_m0, dof, beta)

k = m - 1;
u(m,:) =zeros(1,dof);
rho_m1=zeros(dof,1);

for j=2:k
    a_sum = dot(v(j,:),(f_m - f_m0).');
    u(m,:) = u(m,:)-a_sum*u(j,:);
end
u(m,:)= u(m,:)+beta*(f_m - f_m0).';
u(m,:)= u(m,:)+ rho_m.' - rho_m0.';

for n=2:m
    c_sum = dot(v(n,:),f_m');
    rho_m1 = rho_m1- c_sum*u(n,:).';
end

rho_m1= rho_m1+ beta*f_m;

rho_m1= rho_m1+ rho_m;

rho_m1= abs(rho_m1);

return;