function [Fm,vv,uu,rho,rho0]=new_rho(err_eigvalue, dof, Fm, eigenvec,self_iter_time,self_iter_max,rho,rho0,uu,vv,shell,err_max)

beta=0.15;

Fm0 = Fm;
rho11 = shell*eigenvec.^2;
Fm = rho11 - rho;

if(self_iter_time == 2)
    uu(2,:) = (rho - rho0 + beta*(Fm - Fm0)).';
end

if(self_iter_time >= 2)
    ff =dot(Fm-Fm0, Fm-Fm0);    
    if(abs(ff)<1.e-15)
        fprintf('error: ff == 0 !\n');
    end    
    vv(self_iter_time,:) = (Fm - Fm0).'/ff;   
end

if( err_eigvalue > err_max && self_iter_time < self_iter_max -1)
    
    rho00 = rho0;
    rho0 = rho;
    
    if(self_iter_time >=3 )      
        [Fm,vv,uu,rho]=broyden(self_iter_time,Fm,Fm0,vv,uu,rho0,rho00,dof,beta);
    end
    
    if(self_iter_time<=1)    
        rho =rho+ beta*Fm;
    end
    
    if(self_iter_time==2)
        c =dot(vv(2,:),Fm.');
        rho=rho +beta*Fm-c*uu(2,:).';
    end
    
else
    rho = shell*eigenvec.*eigenvec;
end

rho=abs(rho);

return