function [mat_stif,mat_mass,mat_H] = eig_gstif0(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements)

% gauss=[-0.5773502692,0.5773502692];
% weight= [1.0,1.0];

gauss = [-0.7745966692, 0, 0.7745966692];
weight = [0.5555555555, 0.8888888889, 0.5555555555];

dof = (nx-1)*(ny-1)*(nz-1);

s_size = 64;                        % sparse size
IA     = ones(1, dof * s_size);     % storage indexs for sparse matrix
JA     = ones(1, dof * s_size);
el_mass = zeros(1, dof * s_size);    % for mass matrix
el_stif = zeros(1, dof * s_size);    % for stiff matrix
el_H = zeros(1, dof * s_size); 

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;
hz = (zr-zl)/nz;

%IEL=2;
IEL = 3;
det = hx*hy*hz/8.0;    
node = [[-1,-1,-1];[1,-1,-1];[-1,1,-1];[1,1,-1];[-1,-1,1];[1,-1,1];[-1,1,1];[1,1,1]];

for nn=1:nx*ny*nz
    for i=1:IEL
        for j=1:IEL
            for k=1:IEL
                x = gauss(i);
                y = gauss(j);
                z = gauss(k);
                const = weight(i)*weight(j)*weight(k)*det;
            
                for kk=1:8
                    sf(kk) = 0.125*(1+x*node(kk,1))*(1+y*node(kk,2))*(1+z*node(kk,3));
                    gdsf(kk,1) = 0.25*node(kk,1)*(1+y*node(kk,2))*(1+z*node(kk,3))/hx;
                    gdsf(kk,2) = 0.25*node(kk,2)*(1+x*node(kk,1))*(1+z*node(kk,3))/hy;
                    gdsf(kk,3) = 0.25*node(kk,3)*(1+x*node(kk,1))*(1+y*node(kk,2))/hz;
                end
                
                xx=0.0;
                yy=0.0;
                zz=0.0;
                
                for kk=1:8
                    m = elements(nn,kk);
                    xx = xx + sf(kk)*nodes(m,1);
                    yy = yy + sf(kk)*nodes(m,2);
                    zz = zz + sf(kk)*nodes(m,3);
                end
                
                val=xx^2+yy^2+zz^2;
                val=val/2;
                
               for ni=1:8
                   for nj=1:8
                       mk = elements(nn,ni);
                       mk = nodes(mk,4);
                       ml = elements(nn,nj);
                       ml = nodes(ml,4);
                       if mk~=-1 && ml~=-1
                           index = s_size * (mk-1) + 8 * (ni-1) + nj;
                           IA(index) = mk;
                           JA(index) = ml;
                           el_mass(index) = el_mass(index) + sf(ni)*sf(nj)*const;
                           el_stif(index) = el_stif(index)...
                               +( 0.5*(gdsf(ni,1)*gdsf(nj,1)+gdsf(ni,2)*gdsf(nj,2)+gdsf(ni,3)*gdsf(nj,3))...
                                  + val*sf(ni)*sf(nj) )*const;
                           el_H(index) = el_H(index)...
                               +( gdsf(ni,1)*gdsf(nj,1)+gdsf(ni,2)*gdsf(nj,2)+gdsf(ni,3)*gdsf(nj,3)+sf(ni)*sf(nj))*const;
                       end
                   end
               end
            end
        end
    end
end
mat_stif = sparse(IA, JA, el_stif);
mat_mass = sparse(IA, JA, el_mass);
mat_H = sparse(IA, JA, el_H);


