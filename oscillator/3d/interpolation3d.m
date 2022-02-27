function [U] = interpolation3d(N1,N2,N3,N4,N5,N6,xl,xr,yl,yr,zl,zr,V)                   

hx = (xr-xl)/N1;
hy = (yr-yl)/N2;
hz = (zr-zl)/N3;

Hx = (xr-xl)/N4;
Hy = (yr-yl)/N5;
Hz = (zr-zl)/N6;

u = zeros(N4+1,N5+1,N6+1);

for k=2:N6
    for j=2:N5
        for i=2:N4
            u(i,j,k) = V((N4-1)*(N5-1)*(k-2)+(N4-1)*(j-2)+i-1);
        end
    end
end

uu = zeros(N1+1,N2+1,N3+1);
for k = 2:N3
    for j = 2:N2
        for i = 2:N1
        l = floor((i-1)*hx/Hx)+1;
        m = floor((j-1)*hy/Hy)+1;
        n = floor((k-1)*hz/Hz)+1;
        uu(i,j,k) = (l-(i-1)*hx/Hx)   * (m-(j-1)*hy/Hy)   * (n-(k-1)*hz/Hz)   * u(l,m,n)      +...
                    (1-l+(i-1)*hx/Hx) * (m-(j-1)*hy/Hy)   * (n-(k-1)*hz/Hz)   * u(l+1,m,n)    +...
                    (l-(i-1)*hx/Hx)   * (1-m+(j-1)*hy/Hy) * (n-(k-1)*hz/Hz)   * u(l,m+1,n)    +...
                    (1-l+(i-1)*hx/Hx) * (1-m+(j-1)*hy/Hy) * (n-(k-1)*hz/Hz)   * u(l+1,m+1,n)  +...
                    (l-(i-1)*hx/Hx)   * (m-(j-1)*hy/Hy)   * (1-n+(k-1)*hz/Hz) * u(l,m,n+1)    +...
                    (1-l+(i-1)*hx/Hx) * (m-(j-1)*hy/Hy)   * (1-n+(k-1)*hz/Hz) * u(l+1,m,n+1)  +...
                    (l-(i-1)*hx/Hx)   * (1-m+(j-1)*hy/Hy) * (1-n+(k-1)*hz/Hz) * u(l,m+1,n+1)  +...
                    (1-l+(i-1)*hx/Hx) * (1-m+(j-1)*hy/Hy) * (1-n+(k-1)*hz/Hz) * u(l+1,m+1,n+1);
        end
    end
end

U = zeros((N1-1)*(N2-1)*(N3-1),1);
for k=2:N3
    for j=2:N2
        for i=2:N1
            U((N1-1)*(N2-1)*(k-2)+(N1-1)*(j-2)+i-1, 1) = uu(i,j,k);
        end
    end
end


