function E_total=energy(nx, ny, nz, xl, xr, yl, yr, zl, zr, eigenval, V_Har, rho, nodes, elements)

% gauss = [-0.7745966692, 0, 0.7745966692];
% weight = [0.5555555555, 0.8888888889, 0.5555555555];
gauss  = [-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116]; 
weight = [0.3478548451, 0.6521451548, 0.6521451548, 0.3478548451];

E_total = 0.0;

Ecu = 0;
PE  = 0;
PV  = 0.0;

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;
hz = (zr-zl)/nz;

IEL = 4;
det = hx*hy*hz/8.0;
node = [[-1,-1,-1];[1,-1,-1];[-1,1,-1];[1,1,-1];[-1,-1,1];[1,-1,1];[-1,1,1];[1,1,1]];

for n=1:nx*ny*nz
    for ii=1:8
        imi = elements(n,ii);
        ini = nodes(imi,4);
        if ini~=-1
            el_rho(ii) = rho(ini);
        else
            el_rho(ii) = 0.0;
        end
        
        el_V_H(ii)=V_Har(imi);
    end
    
    Ecu_el = 0;
    PE_el = 0;
    PV_el = 0.0;
    
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
                
                erho = 0.;
                eV_H = 0.;
                for ii=1:8
                    erho = erho +  el_rho(ii) * sf(ii);
                    eV_H = eV_H +  el_V_H(ii) * sf(ii);
                end
                
                vcu_rho_value = eV_H * erho;
                PotExc = -9.0/8*0.77298*(3.0/pi*erho)^(1.0/3.0);
                rho_epsxc_value = erho * PotExc;
                PotVxc = -1.5 * 0.77298 * ((3.0/pi * erho)^(1.0/3.0));
                rho_vxc_value = erho * PotVxc;
                
                Ecu_el = Ecu_el+vcu_rho_value * const;
                PE_el  = PE_el+rho_epsxc_value * const;
                PV_el  = PV_el+rho_vxc_value * const;
            end
        end
    end
    
    E_el = -.5*(Ecu_el) + (PE_el) - (PV_el);
    
    E_total = E_total + E_el;
    Ecu = Ecu + 0.5*Ecu_el;
    PE = PE + PE_el;
    PV = PV + PV_el;
end

E_total=E_total + 2.0*eigenval;

fprintf('Etotal = %f, Ecu = %f, PE = %f, PV = %f\n', E_total, Ecu, PE, PV);

return;