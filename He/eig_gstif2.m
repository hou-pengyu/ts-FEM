function mat_stif = eig_gstif2(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements, V_Har, rho_in)

% gauss = [-0.7745966692, 0, 0.7745966692];
% weight = [0.5555555555, 0.8888888889, 0.5555555555];
gauss  = [-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116]; 
weight = [0.3478548451, 0.6521451548, 0.6521451548, 0.3478548451];

% mat_stif = sparse((nx-1)*(ny-1)*(nz-1), (nx-1)*(ny-1)*(nz-1));
dof = (nx-1)*(ny-1)*(nz-1);

s_size  = 64;                        % sparse size
IA      = ones(1, dof * s_size);     % storage indexs for sparse matrix
JA      = ones(1, dof * s_size);
el_stif = zeros(1, dof * s_size);    % for stiff matrix

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;
hz = (zr-zl)/nz;

IEL = 4;
det = hx*hy*hz/8.0;    
node = [[-1,-1,-1];[1,-1,-1];[-1,1,-1];[1,1,-1];[-1,-1,1];[1,-1,1];[-1,1,1];[1,1,1]];

for nn=1:nx*ny*nz
%     el_stif = zeros(8,8);
    el_rho=zeros(8,1);
    el_V_H=zeros(8,1);
    
    for ii=1:8
        imi = elements(nn,ii);
        ini = nodes(imi,4);
        if ini~=-1
            el_rho(ii) = rho_in(ini);
        else
            el_rho(ii) = 0.0;
        end
        
        el_V_H(ii) = V_Har(imi);
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
                    gdsf(kk,1) = 0.25*node(kk,1)*(1+y*node(kk,2))*(1+z*node(kk,3))/hx;
                    gdsf(kk,2) = 0.25*node(kk,2)*(1+x*node(kk,1))*(1+z*node(kk,3))/hy;
                    gdsf(kk,3) = 0.25*node(kk,3)*(1+x*node(kk,1))*(1+y*node(kk,2))/hz;
                end

                erho = 0.0;
                eV_H = 0.0;
                for kk=1:8
                    erho = erho + sf(kk)*el_rho(kk); %Fang
                    eV_H = eV_H + sf(kk)*el_V_H(kk); %Fang
                end
                eV_xc = -1.5 * 0.77298 * ((3.0/pi * erho)^(1.0/3.0));
                
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
                                +( eV_H + eV_xc )* sf(ni)*sf(nj)*const;       
                         end
                    end
                end
            end
        end
    end
end
mat_stif = sparse(IA, JA, el_stif);


return