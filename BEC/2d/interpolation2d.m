function [U] = interpolation2d(N1,N2,N3,N4,xl,xr,yl,yr,V)                   

hx = (xr-xl)/N1;
hy = (yr-yl)/N2;
Hx = (xr-xl)/N3;
Hy = (yr-yl)/N4;

u = zeros(N3+1,N4+1);

for l=2:N4
    for k=2:N3
        u(k,l)=V((N3-1)*(l-2) + k-1);
    end
end

uu = zeros(N1+1,N2+1);
for j=2:N2
    for i=2:N1
        k = floor((i-1)*hx/Hx)+1;
        l = floor((j-1)*hy/Hy)+1;
        uu(i,j) = (k-(i-1)*hx/Hx)   * (l-(j-1)*hy/Hy)   * u(k,l) +...
                  (1-k+(i-1)*hx/Hx) * (l-(j-1)*hy/Hy)   * u(k+1,l)+...
                  (k-(i-1)*hx/Hx)   * (1-l+(j-1)*hy/Hy) * u(k,l+1)+...
                  (1-k+(i-1)*hx/Hx) * (1-l+(j-1)*hy/Hy) * u(k+1,l+1);
    end
end
U = zeros((N1-1)*(N2-1),1);
for j=2:N2
    for i=2:N1
        U((N1-1)*(j-2) + i-1, 1) = uu(i,j);
    end
end
