function E = energy(nx, ny, xl, xr, yl, yr, nodes, elements, u)

gauss = [-0.7745966692, 0, 0.7745966692];
weight = [0.5555555555, 0.8888888889, 0.5555555555];

beta=200.0;

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;

IEL = 3;
det = hx*hy/4.0;    
node = [[-1,-1];[1,-1];[-1,1];[1,1]];
E = 0.0;
mb = 0.0;
for nn=1:nx*ny
    
    for ii=1:4
        imi = elements(nn,ii);
        ini = nodes(imi,3);
        if ini~=-1
            el_u(ii) = u(ini);
        else
            el_u(ii) = 0.0;
        end
        
    end
    
    E_el = 0;
    mb_el = 0;
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
            
            xx=0.0;
            yy=0.0;
            uu = 0.0;
            grad1 = 0.0;
            grad2 = 0.0;
            for kk=1:4
                m = elements(nn,kk);
                xx = xx + sf(kk)*nodes(m,1);
                yy = yy + sf(kk)*nodes(m,2);
                uu = uu + sf(kk)*el_u(kk);
                grad1 = grad1 + gdsf(kk,1)*el_u(kk);
                grad2 = grad2 + gdsf(kk,2)*el_u(kk);
            end
            
            val1=(xx^2+yy^2)/2.0;
            val2=4.0*exp(-((xx-1.0)^2+yy^2));
            val=val1+val2;
            
            E_el = E_el + (0.5*(grad1^2+grad2^2) + val*uu^2 + 0.5*beta*uu^4)*const;
            mb_el = mb_el + 0.5*beta*uu^4*const;
        end
    end
    E = E + E_el;
    mb = mb + mb_el;
end
mu = E+mb;
fprintf('E = %f, mu = %f\n', E, mu);

return