function x=cg_solve(A, b)

% Solve the linear system Ax=b by the conjugate gradient algorithm

k = 0;
n=size(b,1);
x = zeros(n,1);
p = zeros(n,1);
r = zeros(n,1);
ap = zeros(n,1);

r = A*x;
r = b-r;

while r'*r > 1.0e-20
    k = k+1;
    
    if k==1
        p = r;
    else
        beta = rr1/rr0;
        p = r+beta*p;         
    end
    
    rr0 = r'*r;
    ap = A*p;  
    pap = p'*ap;
    rf = rr0/pap;
    
    x = x+rf*p;
    r = r-rf*ap;
    rr1 = r'*r;
end
  