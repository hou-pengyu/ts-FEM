function [mat_stif, mat_mass] = eig_gstif3(nx, ny, xl, xr, yl, yr, rho_in)

[nodes, elements] = initialize(nx, ny, xl, xr, yl, yr);

gauss = [-0.7745966692, 0, 0.7745966692];
weight = [0.5555555555, 0.8888888889, 0.5555555555];
dof = (nx-1)*(ny-1);

s_size = 16;                        % sparse size
IA     = ones(1, dof * s_size);     % storage indexs for sparse matrix
JA     = ones(1, dof * s_size);
el_mass = zeros(1, dof * s_size);    % for mass matrix
el_stif = zeros(1, dof * s_size);    % for stiff matrix

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;

IEL = 3;
det = hx*hy/4.0;    
node = [[-1,-1];[1,-1];[-1,1];[1,1]];

for nn=1:nx*ny
    
    for ii=1:4
        imi = elements(nn,ii);
        ini = nodes(imi,3);
        if ini~=-1
            el_rho(ii) = rho_in(ini);
        else
            el_rho(ii) = 0.0;
        end
        
    end
    
    for i=1:IEL
        for j=1:IEL
            x = gauss(i);
            y = gauss(j);
            const = weight(i)*weight(j)*det;
            
            for kk=1:4
                sf(kk) = 0.25*(1+x*node(kk,1))*(1+y*node(kk,2));
                gdsf(kk,1) = 0.5*node(kk,1)*(1+y*node(kk,2))/hx;
                gdsf(kk,2) = 0.5*node(kk,2)*(1+x*node(kk,1))/hy;

            end
            
            erho = 0.0;
            xx=0.0;
            yy=0.0;            
            for kk=1:4
                m = elements(nn,kk);
                xx = xx + sf(kk)*nodes(m,1);
                yy = yy + sf(kk)*nodes(m,2);   
                erho = erho + sf(kk)*el_rho(kk);
            end
            val=xx^2+yy^2;
            val=val/2;
            
            for ni=1:4
                for nj=1:4
                    mk = elements(nn,ni);
                    mk = nodes(mk,3);
                    ml = elements(nn,nj);
                    ml = nodes(ml,3);
                    if mk~=-1 && ml~=-1
                        index = s_size * (mk-1) + 4 * (ni-1) + nj;
                        IA(index) = mk;
                        JA(index) = ml;
                        el_mass(index) = el_mass(index) + sf(ni)*sf(nj)*const;
                        el_stif(index) = el_stif(index)...
                             +( 0.5*(gdsf(ni,1)*gdsf(nj,1)+gdsf(ni,2)*gdsf(nj,2))...
                            + (val+erho)*sf(ni)*sf(nj) )*const;                          
                    end
                end
            end
            
        end
    end
end
mat_stif = sparse(IA, JA, el_stif);
mat_mass = sparse(IA, JA, el_mass);

return