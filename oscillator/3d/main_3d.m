% Harmonic oscillator
% 3d nonlinear eigenvalue problems

xl = -5;
xr = -xl;
yl = xl;
yr = xr;
zl = xl;
zr = xr; % the bound of the computational domain

nx = 160;
ny = nx;
nz = nx; % mesh of reference solution
t_lambda = 0.0;
% [A,B,H,D,U] = impl_oscillator(nx, ny, nz, xl, xr, yl, yr, zl, zr, t_lambda);  % reference solution
load('B_160_5.mat'); % mass matrix
load('D_160_5.mat'); % reference eigenvalue
load('H_160_5.mat');
load('U_160_5.mat'); % reference eigenvector

N1 = [32,40,48,58,68,78];  % mesh
N2 = [18,20,22,24,26,28];  

l = length(N1);
h  = zeros(l,1);

l_h      = zeros(l,1);
l_x      = zeros(l,1);
l_y      = zeros(l,1);
l_z      = zeros(l,1);
l_H      = zeros(l,1);
l_ts     = zeros(l,1);
l_ts_lin = zeros(l,1);

l_err_fem    = zeros(l,1);
l_err_ts     = zeros(l,1);
l_err_ts_lin = zeros(l,1);

H1_err_fem    = zeros(l,1); 
H1_err_ts     = zeros(l,1);
H1_err_ts_lin = zeros(l,1);

time_hhh    = zeros(l,1);
time_fem    = zeros(l,1);

time_HHH    = zeros(l,1);
time_hHH    = zeros(l,1);
time_HhH    = zeros(l,1);
time_HHh    = zeros(l,1);
time_uHHH_h = zeros(l,1);
time_lamHHH_h  = zeros(l,1);
time_ts     = zeros(l,1);

time_ini_lin   = zeros(l,1);
time_hHH_lin   = zeros(l,1);
time_HhH_lin   = zeros(l,1);
time_HHh_lin   = zeros(l,1);
time_uHHH_h_lin = zeros(l,1);
time_lamHHH_h_lin = zeros(l,1);
time_ts_lin    = zeros(l,1);

for i=1:l
    fprintf('%i-th\n', i)
    fprintf('FEM\n')
    h(i) = (xr-xl)/N1(i); % step size
    
    %========================================================================%
    % solve nonlinear eigenvalue problems by classic FEM

    tic;
    [mat_stif_h, mat_mass_h, mat_H_h, l_h(i,1), u_h] = impl_oscillator(N1(i),N1(i),N1(i), xl, xr, yl, yr, zl, zr, t_lambda); % solutions of FEM
    time_fem(i,1) = toc;
    
    tic
    uhhh = interpolation3d(nx, ny, nz, N1(i), N1(i), N1(i), xl, xr, yl, yr, zl, zr, u_h);
    if uhhh'*U<0
        uhhh = -uhhh;
    end
    time_hhh(i,1) = toc;
        
    uhhh_err  = uhhh - U;
    H1_err_fem(i) = sqrt(uhhh_err'*H*uhhh_err); % error of eigenfunction of FEM in H1 norm
    l_err_fem(i,1) = abs(l_h(i,1) - D);         % error of eigenvalue of FEM

    %========================================================================%
    % solve nonlinear eigenvalue problems by basic two-scale FEM
    fprintf('%i-th\n', i)
    fprintf('two-scale FEM\n')    
    
    tic % solutions of l_HHH and u_HHH
    [mat_stif_H, mat_mass_H, mat_H_H, l_H(i,1), u_H] = impl_oscillator(N2(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, t_lambda);
    time_HHH(i,1) = toc;
    
    tic % solutions of l_hHH and u_hHH
    [mat_stif_x, mat_mass_x, mat_H_x, l_x(i,1), u_x] = impl_oscillator(N1(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, t_lambda);
    time_hHH(i,1) = toc;
    
    tic % solutions of l_HhH and u_HhH
    [mat_stif_y, mat_mass_y, mat_H_y, l_y(i,1), u_y] = impl_oscillator(N2(i), N1(i), N2(i), xl, xr, yl, yr, zl, zr, t_lambda);
    time_HhH(i,1) = toc;
    
    tic % solutions of l_HHh and u_HHh
    [mat_stif_z, mat_mass_z, mat_H_z, l_z(i,1), u_z] = impl_oscillator(N2(i), N2(i), N1(i), xl, xr, yl, yr, zl, zr, t_lambda);
    time_HHh(i,1) = toc;
    
    tic
    uHHH = interpolation3d(N1(i), N1(i), N1(i), N2(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, u_H);
    uhHH = interpolation3d(N1(i), N1(i), N1(i), N1(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, u_x);
    uHhH = interpolation3d(N1(i), N1(i), N1(i), N2(i), N1(i), N2(i), xl, xr, yl, yr, zl, zr, u_y);
    uHHh = interpolation3d(N1(i), N1(i), N1(i), N2(i), N2(i), N1(i), xl, xr, yl, yr, zl, zr, u_z);
    
%     if uHHH'*u_h<0
%         uHHH = -uHHH;
%     end
%     if uhHH'*u_h<0
%         uhHH = -uhHH;
%     end
%     if uHhH'*u_h<0
%         uHhH = -uHhH;
%     end
%     if uHHh'*u_h<0
%         uHHh = -uHHh;
%     end

    if uhHH'*uHHH<0
        uhHH = -uhHH;
    end
    if uHhH'*uHHH<0
        uHhH = -uHhH;
    end
    if uHHh'*uHHH<0
        uHHh = -uHHh;
    end
    
    u_ts = uhHH + uHhH + uHHh - 2*uHHH;  % eigenvector of two-scale FEM
    mod_u_ts = module_u(N1(i), N1(i), N1(i), xl, xr, yl, yr, zl, zr,u_ts);
    u_ts = u_ts/mod_u_ts;

    time_uHHH_h(i,1) = toc;
            
    tic
    [mat_stif_u_ts, mat_mass_u_ts] = eig_gstif3(N1(i), N1(i), N1(i), xl, xr, yl, yr, zl, zr, u_ts.^2);
    l_ts(i,1) = (u_ts'*(mat_stif_u_ts)*u_ts)/(u_ts'*mat_mass_u_ts*u_ts); % eigenvalue of two-scale FEM
    time_lamHHH_h(i,1) = toc;
        
    time_ts(i,1)  = time_HHH(i,1) + time_hHH(i,1) + time_HhH(i,1) + time_HHh(i,1) +...
                    time_uHHH_h(i,1) + time_lamHHH_h(i,1) ; % time of two-scale FEM 
                
    u_ts_h = interpolation3d(nx, ny, nz, N1(i), N1(i), N1(i), xl, xr, yl, yr, zl, zr, u_ts);
    if u_ts_h'*U<0
        u_ts_h = -u_ts_h;
    end
    u_ts_h_err = u_ts_h - U;
    H1_err_ts(i) = sqrt(u_ts_h_err'*H*u_ts_h_err); % error of eigenfunction of two-scale FEM in H1 norm
    l_err_ts(i,1)  = abs(l_ts(i,1) - D);           % error of eigenvalue of two-scale FEM

    %========================================================================%
    % solve nonlinear eigenvalue problems by linearized two-scale FEM
    fprintf('%i-th\n', i)
    fprintf('linearized two-scale FEM\n')  
    
    tic    
    u_hHH = interpolation3d(N1(i), N2(i), N2(i), N2(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, u_H);
    u_HhH = interpolation3d(N2(i), N1(i), N2(i), N2(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, u_H);
    u_HHh = interpolation3d(N2(i), N2(i), N1(i), N2(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, u_H);
    time_ini_lin(i,1) = toc;
    
    tic % solution of u^hHH
    u_x_lin = impl_oscillator_linearized(N1(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, u_hHH, l_H(i,1), t_lambda);
    time_hHH_lin(i,1) = toc;
    
    tic % solution of u^HhH
    u_y_lin = impl_oscillator_linearized(N2(i), N1(i), N2(i), xl, xr, yl, yr, zl, zr, u_HhH, l_H(i,1), t_lambda);
    time_HhH_lin(i,1) = toc;
    
    tic % solution of u^HHh
    u_z_lin = impl_oscillator_linearized(N2(i), N2(i), N1(i), xl, xr, yl, yr, zl, zr, u_HHh, l_H(i,1), t_lambda);
    time_HHh_lin(i,1) = toc;
    
    tic
    uhHH_lin = interpolation3d(N1(i), N1(i), N1(i), N1(i), N2(i), N2(i), xl, xr, yl, yr, zl, zr, u_x_lin);
    uHhH_lin = interpolation3d(N1(i), N1(i), N1(i), N2(i), N1(i), N2(i), xl, xr, yl, yr, zl, zr, u_y_lin);
    uHHh_lin = interpolation3d(N1(i), N1(i), N1(i), N2(i), N2(i), N1(i), xl, xr, yl, yr, zl, zr, u_z_lin);
    
%     if uhHH_lin'*u_h<0
%         uhHH_lin = -uhHH_lin;
%     end
%     if uHhH_lin'*u_h<0
%         uHhH_lin = -uHhH_lin;
%     end
%     if uHHh_lin'*u_h<0
%         uHHh_lin = -uHHh_lin;
%     end

    if uhHH_lin'*uHHH<0
        uhHH_lin = -uhHH_lin;
    end
    if uHhH_lin'*uHHH<0
        uHhH_lin = -uHhH_lin;
    end
    if uHHh_lin'*uHHH<0
        uHHh_lin = -uHHh_lin;
    end
    
    u_ts_lin = uhHH_lin + uHhH_lin + uHHh_lin - 2*uHHH; % eigenvector of linearized two-scale FEM
    mod_u_ts_lin = module_u(N1(i), N1(i), N1(i), xl, xr, yl, yr, zl, zr, u_ts_lin);
    u_ts_lin = u_ts_lin/mod_u_ts_lin;

    time_uHHH_h_lin(i,1) = toc;
    
    tic
    [mat_stif_u_ts_lin, mat_mass_u_ts_lin] = eig_gstif3(N1(i), N1(i), N1(i), xl, xr, yl, yr, zl, zr, u_ts_lin.^2);
    l_ts_lin(i,1)  = (u_ts_lin'*(mat_stif_u_ts_lin)*u_ts_lin)/(u_ts_lin'*mat_mass_u_ts_lin*u_ts_lin); % eigenvalue of linearized two-scale FEM
    time_lamHHH_h_lin(i,1) = toc;
       
    time_ts_lin(i,1) = time_HHH(i,1) + time_ini_lin(i,1) + time_hHH_lin(i,1) + time_HhH_lin(i,1)...
                      + time_HHh_lin(i,1) + time_uHHH_h_lin(i,1) + time_lamHHH_h_lin(i,1); % time of inearized two-scale FEM
                  
    u_ts_lin_h = interpolation3d(nx, ny, nz, N1(i), N1(i), N1(i), xl, xr, yl, yr, zl, zr, u_ts_lin);
    if u_ts_lin_h'*U<0
        u_ts_lin_h = -u_ts_lin_h;
    end
    u_ts_lin_h_err = u_ts_lin_h - U;
    H1_err_ts_lin(i) = sqrt(u_ts_lin_h_err' * H * u_ts_lin_h_err); % error of eigenfunction of linearized two-scale FEM in H1 norm  
    l_err_ts_lin(i,1) = abs(l_ts_lin(i,1) - D);                    % error of eigenvalue of linearized two-scale FEM

end

figure(1) % The convergence curves of eigenvalue
loglog(h,l_err_fem(:,1),'mo-','LineWidth', 1.8,'MarkerSize',15);
hold on
loglog(h,l_err_ts(:,1),'c^--','LineWidth', 1.8,'MarkerSize',12);
hold on
loglog(h,l_err_ts_lin(:,1),'ks-.','LineWidth', 1.8,'MarkerSize',7);
hold on

x=zeros();
y=x;
for i=1:100
    x(i)=0.01*i;
    y(i)=x(i).^2/10;

end
loglog(x,y,'b-.','LineWidth', 1.8);
hold on

xlabel({'$h$'},'Interpreter','latex');
ylabel('error');
legend({'$|\lambda_{ref} - \lambda_{h,h,h}|$',...
        '$|\lambda_{ref} - \lambda_{H,H,H}^{h}|$',...
        '$|\lambda_{ref} - \widetilde{\lambda}_{H,H,H}^{h}|$','slope=2'},'Interpreter','latex');
set(gca,'XLim',[0.12,max(h)+0.02]); 
set(gca,'YLim',[7*10^(-4),1.1*10^(-2)]);
set(gca,'FontSize',25);

figure(2) % The convergence curves of eigenfunction
loglog(h,H1_err_fem,'mo-','linewidth',1.8,'MarkerSize',15);
hold on
loglog(h,H1_err_ts,'c^--','linewidth',1.8,'MarkerSize',12);
hold on
loglog(h,H1_err_ts_lin,'ks-.','linewidth',1.8,'MarkerSize',7);
hold on

x=zeros();
y=x;
for i=1:100
    x(i)=0.01*i;
    y(i)=x(i)/2.2;
end
loglog(x,y,'b-.','linewidth',1.8);
hold on

xlabel({'$h$'},'Interpreter','latex');
ylabel('error');
legend({'$\| u_{ref} - u_{h,h,h}\|_{1}$',...
        '$\| u_{ref} - u_{H,H,H}^{h}\|_{1}$',...
        '$\| u_{ref} - \widetilde{u}_{H,H,H}^{h}\|_{1}$ ','slope=1'},'Interpreter','latex');
set(gca,'XLim',[0.12,max(h)+0.02]);    
set(gca,'YLim',[3*10^(-2),1.5*10^(-1)]);

set(gca,'YTick',10^(-1));
set(gca,'FontSize',25);

figure(3) % The running time
dh = (xr-xl)./h;
plot(dh,time_fem(:,1),'mo-','Linewidth',1.8,'MarkerSize',15);
hold on
plot(dh,time_ts(:,1),'c^--','Linewidth',1.8,'MarkerSize',12);
hold on
plot(dh,time_ts_lin(:,1),'ks-.','Linewidth',1.8,'MarkerSize',7);
hold on

xlabel({'$L/h$'},'Interpreter','latex');
ylabel('time(s)');
legend('Standard FEM','Algorithm 5.1','Algorithm 5.2');
set(gca,'XLim',[30.0,80.0]);   
set(gca,'YLim',[-200,time_fem(l,1)+300]);

set(gca,'FontSize',25);

