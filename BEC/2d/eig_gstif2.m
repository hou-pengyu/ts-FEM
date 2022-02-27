function mat_stif = eig_gstif2(nx, ny, xl, xr, yl, yr, nodes, elements, rho_in)

gauss = [-0.7745966692, 0, 0.7745966692];
weight = [0.5555555555, 0.8888888889, 0.5555555555];
dof = (nx-1)*(ny-1);

s_size = 16;                        % sparse size
IA     = ones(1, dof * s_size);     % storage indexs for sparse matrix
JA     = ones(1, dof * s_size);
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
            end
            
            erho = 0.0;
            
            for kk=1:4
                
                erho = erho + sf(kk)*el_rho(kk);
                
            end
            
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
                        el_stif(index) = el_stif(index)...
                            +(200*erho)* sf(ni)*sf(nj)*const;
                    end
                end
            end
            
        end
    end
end
mat_stif = sparse(IA, JA, el_stif);

return