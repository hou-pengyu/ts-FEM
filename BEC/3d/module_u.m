function mod_u = module_u(nx, ny, nz, xl, xr, yl, yr, zl, zr, u)

gauss = [-0.7745966692, 0, 0.7745966692];
weight = [0.5555555555, 0.8888888889, 0.5555555555];

nodes = zeros((nx+1)*(ny+1)*(nz+1), 4);
elements = zeros(nx*ny*nz, 8);

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;
hz = (zr-zl)/nz;
for i=1:nx+1
    for j=1:ny+1
        for k=1:nz+1
            m = (nx+1)*(ny+1)*(k-1) + (nx+1)*(j-1) + i;
            nodes(m,1) = xl + hx*(i-1);
            nodes(m,2) = yl + hy*(j-1);
            nodes(m,3) = zl + hz*(k-1);
            nodes(m,4) = m;
        end
    end
end

for i=1:nx
    for j=1:ny
        for k=1:nz
            m = nx*ny*(k-1) + nx*(j-1) + i;
            elements(m,1) = (nx+1)*(ny+1)*(k-1)+(nx+1)*(j-1)+i;
            elements(m,2) = elements(m,1)+1;
            elements(m,3) = elements(m,1)+nx+1;
            elements(m,4) = elements(m,3)+1;  
            elements(m,5) = elements(m,1)+(nx+1)*(ny+1);
            elements(m,6) = elements(m,5)+1;
            elements(m,7) = elements(m,5)+nx+1;
            elements(m,8) = elements(m,7)+1;
        end
    end
end

k = 0;
for i=1:(nx+1)*(ny+1)*(nz+1)
    bound_val = (nodes(i,1)-xl)*(nodes(i,1)-xr)*(nodes(i,2)-yl)*(nodes(i,2)-yr)*(nodes(i,3)-zl)*(nodes(i,3)-zr);
    if bound_val==0
        k = k+1;
        nodes(i,4) = -1;
    else
        nodes(i,4) = nodes(i,4)-k;
    end    
end

IEL = 3;
det = hx*hy*hz/8.0;    
node = [[-1,-1,-1];[1,-1,-1];[-1,1,-1];[1,1,-1];[-1,-1,1];[1,-1,1];[-1,1,1];[1,1,1]];

mod_u = 0.0;
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
                end
                
                u0=0.0;
                for kk=1:8
                    m = elements(nn,kk);
                    m = nodes(m,4);
                    if m~=-1
                        u0 = u0 + u(m)*sf(kk);
                    end
                end
                
                mod_u = mod_u + u0^2 * const;
            end
        end
    end
end

mod_u = sqrt(mod_u);
 