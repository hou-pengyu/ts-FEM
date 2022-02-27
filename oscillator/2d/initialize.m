function [nodes, elements] = initialize(nx, ny, xl, xr, yl, yr)

gauss = [-0.7745966692, 0, 0.7745966692];
weight = [0.5555555555, 0.8888888889, 0.5555555555];

nodes = zeros((nx+1)*(ny+1), 3);
elements = zeros(nx*ny, 4);

hx = (xr-xl)/nx;
hy = (yr-yl)/ny;
for i=1:nx+1
    for j=1:ny+1
        m = (nx+1)*(j-1) + i;
        nodes(m,1) = xl + hx*(i-1);
        nodes(m,2) = yl + hy*(j-1);
        nodes(m,3) = m;
    end
end

for i=1:nx
    for j=1:ny
        m = nx*(j-1) + i;
        elements(m,1) = (nx+1)*(j-1)+i;
        elements(m,2) = elements(m,1)+1;
        elements(m,3) = elements(m,1)+nx+1;
        elements(m,4) = elements(m,3)+1;
    end
end

k = 0;
for i=1:(nx+1)*(ny+1)
    bound_val = (nodes(i,1)-xl)*(nodes(i,1)-xr)*(nodes(i,2)-yl)*(nodes(i,2)-yr);
    if bound_val==0
        k = k+1;
        nodes(i,3) = -1;
    else
        nodes(i,3) = nodes(i,3)-k;
    end    
end
return