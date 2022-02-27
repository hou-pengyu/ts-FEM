function [V_Har] = V_Har_pot(nx, ny, nz, xl, xr, yl, yr, zl, zr, nodes, elements, rho_in, mat_stif)

% gauss = [-0.7745966692, 0, 0.7745966692];
% weight = [0.5555555555, 0.8888888889, 0.5555555555];
gauss  = [-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116]; 
weight = [0.3478548451, 0.6521451548, 0.6521451548, 0.3478548451];

% f = zeros((nx-1)*(ny-1)*(nz-1),1);
dof = (nx-1)*(ny-1)*(nz-1);
IA     = ones(1, dof * 8);     % storage indexs 
el_f   = zeros(1, dof * 8);    % for f

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;
hz = (zr-zl)/nz;

IEL = 4;
det = hx*hy*hz/8.0;    
node = [[-1,-1,-1];[1,-1,-1];[-1,1,-1];[1,1,-1];[-1,-1,1];[1,-1,1];[-1,1,1];[1,1,1]];

for nn=1:nx*ny*nz
%     el_f = zeros(8,1);
    el_rho=zeros(8,1);
    for ii=1:8
        imi = elements(nn,ii); %all of nodes including nodes on the boundary
        ini = nodes(imi,4);    %indicating whether a point is on the boundary or not
        if(ini > 0)
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
                
                xx=0.0;
                yy=0.0;
                zz=0.0;
                rho = 0.0;
                for kk=1:8
                    m = elements(nn,kk);
                    xx = xx + sf(kk)*nodes(m,1);
                    yy = yy + sf(kk)*nodes(m,2);
                    zz = zz + sf(kk)*nodes(m,3);
                    
                    rho = rho + sf(kk)*el_rho(kk);
                end
                rr = sqrt(xx^2+yy^2+zz^2);
%                 
                for ni = 1:8
                    mk = elements(nn,ni);
                    mk = nodes(mk,4);
                    if mk~=-1
                       index = 8*(mk-1) + ni;
                       IA(index) = mk;
                       if rr>=4
                          el_f(index) = el_f(index) + 4*pi*rho*sf(ni) * const;
                       else
                          el_f(index) = el_f(index) + (4*pi*rho-3.0/32.0)*sf(ni) * const;
                       end
                    end
                end
%                 for ni=1:8
%                     if rr>=4 
%                         el_f(ni) = el_f(ni) + 4*pi*rho*sf(ni) * const;
%                     else 
%                         el_f(ni) = el_f(ni) + (4*pi*rho-3.0/32.0)*sf(ni) * const;
%                     end
%                 end
            end
        end
    end
%     for ni=1:8
%         mi = elements(nn,ni);
%         mi = nodes(mi,4);
%         if mi~=-1
%             f(mi) = f(mi)+el_f(ni);            
%         end
%     end
end
f = sparse(IA, 1, el_f);
mat_stif=2*mat_stif;
uu = cg_solve(mat_stif, f);

V_Har = zeros((nx+1)*(ny+1)*(nz+1), 1);
RR = 4.0;
for i=1:nx+1
    for j=1:ny+1
        for k=1:nz+1
            m = (nx+1)*(ny+1)*(k-1) + (nx+1)*(j-1) + i;
            n = nodes(m,4);
            x = nodes(m,1);
            y = nodes(m,2);
            z = nodes(m,3);
            val = sqrt(x^2+y^2+z^2);

            if n == -1
                if val>RR
                   V_Har(m) = 2/val;                
                else
                   V_Har(m) = 2*(3*RR^2-val^2)/(2*RR^3);
                end
            else 
                if val>RR
                   V_Har(m) = uu(n)+2/val;
                else
                   V_Har(m) = uu(n)+2*(3*RR^2-val^2)/(2*RR^3);
                end
            end
        end
    end
end
return

