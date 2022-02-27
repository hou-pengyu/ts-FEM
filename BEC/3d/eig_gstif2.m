function mat_stif = eig_gstif2(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements, rho_in)

gauss = [-0.7745966692, 0, 0.7745966692];
weight = [0.5555555555, 0.8888888889, 0.5555555555];

% mat_stif = sparse((nx-1)*(ny-1)*(nz-1), (nx-1)*(ny-1)*(nz-1));
dof = (nx-1)*(ny-1)*(nz-1);

s_size = 64;                        % sparse size
IA     = ones(1, dof * s_size);     % storage indexs for sparse matrix
JA     = ones(1, dof * s_size);
el_stif = zeros(1, dof * s_size);    % for stiff matrix

beta=200.0;

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;
hz = (zr-zl)/nz;

IEL = 3;
det = hx*hy*hz/8.0;    
node = [[-1,-1,-1];[1,-1,-1];[-1,1,-1];[1,1,-1];[-1,-1,1];[1,-1,1];[-1,1,1];[1,1,1]];

for nn=1:nx*ny*nz
%     el_stif = zeros(8,8);
%     el_rho=zeros(8,1);
   
    for ii=1:8
        imi = elements(nn,ii);
        ini = nodes(imi,4);
        if ini~=-1
            el_rho(ii) = rho_in(ini);
        else
            el_rho(ii) = 0.0;
        end
        
    end
       
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
                
                erho = 0.0;               
                for kk=1:8           
                    erho = erho + sf(kk)*el_rho(kk); 
                end
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
                            el_stif(index) = el_stif(index)...
                                +( beta*erho )* sf(ni)*sf(nj)*const;          
                         end
                    end
                end
%                 for ni=1:8
%                     for nj=1:8
%                         el_stif(ni,nj) = el_stif(ni,nj)...
%                             + (beta*erho) * sf(ni)*sf(nj)*const; 
%                     end
%                 end
            end
        end
    end
%     for ni=1:8
%         for nj=1:8
%             mi = elements(nn,ni);
%             mi = nodes(mi,4);
%             mj = elements(nn,nj);
%             mj = nodes(mj,4);
%             if mi~=-1 && mj~=-1
%                 mat_stif(mi,mj) = mat_stif(mi,mj) + el_stif(ni,nj);
%             end
%         end
%     end
end
mat_stif = sparse(IA, JA, el_stif);

return