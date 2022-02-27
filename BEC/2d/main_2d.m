% BEC
% 2d nonlinear eigenvalue problems

xl = -5;
xr = -xl;
yl = xl;
yr = xr; % the bound of the computational domain

nx = 160;
ny = nx; % mesh of reference solution
t_lambda = 9.0;
[~,B,H,D,U] = impl_GPE(nx, ny, xl, xr, yl, yr, t_lambda);  % reference solution

N1 = [32,40,48,58,68,78];   % mesh
N2 = [18,20,22,24,26,28]; 

l = length(N1);
h  = zeros(l,1);

l_h      = zeros(l,1);
l_x      = zeros(l,1);
l_y      = zeros(l,1);
l_H      = zeros(l,1);
l_ts     = zeros(l,1);
l_ts_lin = zeros(l,1);

l_err_fem    = zeros(l,1);
l_err_ts     = zeros(l,1);
l_err_ts_lin = zeros(l,1);

H1_err_fem    = zeros(l,1); 
H1_err_ts     = zeros(l,1);
H1_err_ts_lin = zeros(l,1);

time_hh    = zeros(l,1);
time_fem    = zeros(l,1);

time_HH     = zeros(l,1);
time_hH     = zeros(l,1);
time_Hh     = zeros(l,1);
time_ts     = zeros(l,1);
time_uHH_h  = zeros(l,1);
time_lamHH_h  = zeros(l,1);

time_ini_lin   = zeros(l,1);
time_hH_lin    = zeros(l,1);
time_Hh_lin    = zeros(l,1);
time_ts_lin    = zeros(l,1);
time_uHH_h_lin = zeros(l,1);
time_lamHH_h_lin = zeros(l,1);

for i=1:l
    
    h(i) = (xr-xl)/N1(i); % step size
    %========================================================================%
    % solve nonlinear eigenvalue problems by classic FEM
    fprintf('%i-th\n', i)
    fprintf('FEM\n')
    
    tic;
    [mat_stif_h, mat_mass_h, mat_H_h, l_h(i,1), u_h] = impl_GPE(N1(i), N1(i), xl, xr, yl, yr, t_lambda); % solutions of FEM
    time_fem(i,1) = toc;
    
    tic
    uhh = interpolation2d(nx, ny, N1(i), N1(i), xl, xr, yl, yr, u_h);
    if uhh'*U<0
        uhh = -uhh;
    end
    time_hh(i,1) = toc;
    
    uhh_err  = uhh - U;
    H1_err_fem(i) = sqrt(uhh_err'*H*uhh_err); % error of eigenfunction of FEM in H1 norm
    l_err_fem(i,1) = abs(l_h(i,1) - D);       % error of eigenvalue of FEM

    %===============================================================================================%
    % solve nonlinear eigenvalue problems by basic two-scale FEM
    fprintf('%i-th\n', i)
    fprintf('two-scale FEM\n') 
    
    tic % solutions of l_HH and u_HH
    [mat_stif_H, mat_mass_H, mat_H_H, l_H(i,1), u_H] = impl_GPE(N2(i), N2(i), xl, xr, yl, yr, t_lambda);
    time_HH(i,1) = toc;
    
    tic % solutions of l_hH and u_hH
    [mat_stif_x, mat_mass_x, mat_H_x, l_x(i,1), u_x] = impl_GPE(N1(i), N2(i), xl, xr, yl, yr, t_lambda);
    time_hH(i,1) = toc;
    
    tic % solutions of l_Hh and u_Hh
    [mat_stif_y, mat_mass_y, mat_H_y, l_y(i,1), u_y] = impl_GPE(N2(i), N1(i), xl, xr, yl, yr, t_lambda);
    time_Hh(i,1) = toc;
    
    tic
    uHH = interpolation2d(N1(i), N1(i), N2(i), N2(i), xl, xr, yl, yr, u_H);
    uhH = interpolation2d(N1(i), N1(i), N1(i), N2(i), xl, xr, yl, yr, u_x);
    uHh = interpolation2d(N1(i), N1(i), N2(i), N1(i), xl, xr, yl, yr, u_y);
    
    if uHH'*u_h<0
        uHH = -uHH;
    end
    if uhH'*u_h<0
        uhH = -uhH;
    end
    if uHh'*u_h<0
        uHh = -uHh;
    end
    
    u_ts = uhH + uHh - uHH; % eigenvector of two-scale FEM
    mod_u_ts = module_u(N1(i), N1(i), xl, xr, yl, yr, u_ts);
    u_ts = u_ts/mod_u_ts;

    time_uHH_h(i,1) = toc;
    
    tic
    [mat_stif_u_ts, mat_mass_u_ts] = eig_gstif3(N1(i), N1(i), xl, xr, yl, yr, u_ts.^2);
    l_ts(i,1) = (u_ts'*(mat_stif_u_ts)*u_ts)/(u_ts'*mat_mass_u_ts*u_ts); % eigenvalue of two-scale FEM
    time_lamHH_h(i,1) = toc;
        
    time_ts(i,1)  = time_HH(i,1) + time_hH(i,1) + time_Hh(i,1) +...
                    time_uHH_h(i,1) + time_lamHH_h(i,1); % time of two-scale FEM
                
    u_ts_h = interpolation2d(nx, ny, N1(i), N1(i), xl, xr, yl, yr, u_ts);
    if u_ts_h'*U<0
        u_ts_h = -u_ts_h;
    end
    u_ts_h_err = u_ts_h - U;
    H1_err_ts(i) = sqrt(u_ts_h_err'*H*u_ts_h_err);  % error of eigenfunction of two-scale FEM in H1 norm
    l_err_ts(i,1)  = abs(l_ts(i,1) - D);            % error of eigenvalue of two-scale FEM

    %===============================================================================================%
    % solve nonlinear eigenvalue problems by linearized two-scale FEM
    fprintf('%i-th\n', i)
    fprintf('linearized two-scale FEM\n') 
    
    tic    
    u_hH = interpolation2d(N1(i), N2(i), N2(i), N2(i), xl, xr, yl, yr, u_H);
    u_Hh = interpolation2d(N2(i), N1(i), N2(i), N2(i), xl, xr, yl, yr, u_H);
    time_ini_lin(i,1) = toc;
    
    tic % solution of u^hH
    u_x_lin = impl_GPE_linearized(N1(i), N2(i), xl, xr, yl, yr, u_hH, l_H(i,1), t_lambda);
    time_hH_lin(i,1) = toc;
    
    tic % solution of u^Hh
    u_y_lin = impl_GPE_linearized(N2(i), N1(i), xl, xr, yl, yr, u_Hh, l_H(i,1), t_lambda);
    time_Hh_lin(i,1) = toc;
    
    tic
    uhH_lin = interpolation2d(N1(i), N1(i), N1(i), N2(i), xl, xr, yl, yr, u_x_lin);
    uHh_lin = interpolation2d(N1(i), N1(i), N2(i), N1(i), xl, xr, yl, yr, u_y_lin);

    if uhH_lin'*u_h<0
        uhH_lin = -uhH_lin;
    end
    if uHh_lin'*u_h<0
        uHh_lin = -uHh_lin;
    end
    
    u_ts_lin = uhH_lin + uHh_lin - uHH; % eigenvector of linearized two-scale FEM
    mod_u_ts_lin = module_u(N1(i), N1(i), xl, xr, yl, yr, u_ts_lin);
    u_ts_lin = u_ts_lin/mod_u_ts_lin;

    time_uHH_h_lin(i,1) = toc;
    
    tic
    [mat_stif_u_ts_lin,mat_mass_u_ts_lin] = eig_gstif3(N1(i), N1(i), xl, xr, yl, yr, u_ts_lin.^2);
    l_ts_lin(i,1)  = (u_ts_lin'*(mat_stif_u_ts_lin)*u_ts_lin)/(u_ts_lin'*mat_mass_u_ts_lin*u_ts_lin); % eigenvalue of linearized two-scale FEM
    time_lamHH_h_lin(i,1) = toc;
    
    time_ts_lin(i,1) = time_HH(i,1) + time_ini_lin(i,1) + time_hH_lin(i,1) + time_Hh_lin(i,1)...
                       + time_uHH_h_lin(i,1) + time_lamHH_h_lin(i,1); % time of inearized two-scale FEM
                   
    u_ts_lin_h = interpolation2d(nx, ny, N1(i), N1(i), xl, xr, yl, yr, u_ts_lin);
    if u_ts_lin_h'*U<0
        u_ts_lin_h = -u_ts_lin_h;
    end                  
    u_ts_lin_h_err = u_ts_lin_h - U;
    H1_err_ts_lin(i) = sqrt(u_ts_lin_h_err' * H * u_ts_lin_h_err);  % error of eigenfunction of linearized two-scale FEM in H1 norm
    l_err_ts_lin(i,1) = abs(l_ts_lin(i,1) - D);                     % error of eigenvalue of linearized two-scale FEM 

end

figure(1) % The convergence curves of eigenvalue
loglog(h,l_err_fem(:,1),'mo-','LineWidth', 1.8, 'MarkerSize',15);
hold on
loglog(h,l_err_ts(:,1),'c^--','LineWidth', 1.8, 'MarkerSize',12);
hold on
loglog(h,l_err_ts_lin(:,1),'ks-.','LineWidth', 1.8, 'MarkerSize',7);
hold on

x=zeros();
y=x;
for i=1:100
    x(i)=0.01*i;
    y(i)=x(i).^2/7;

end
loglog(x,y,'b-.','LineWidth', 1.8);
hold on

xlabel({'$h$'},'Interpreter','latex');
ylabel('error');
legend({'$|\lambda_{ref} - \lambda_{h,h}|$',...
        '$|\lambda_{ref} - \lambda_{H,H}^{h}|$',...
        '$|\lambda_{ref} - \widetilde{\lambda}_{H,H}^{h}|$','slope=2'},'Interpreter','latex');
set(gca,'XLim',[0.12,max(h)+0.02]);
set(gca,'YLim',[1.2*10^(-3),2*10^(-2)]);
set(gca,'FontSize',25);

figure(2) % The convergence curves of eigenfunction
loglog(h,H1_err_fem,'mo-','linewidth',1.8, 'MarkerSize',15);
hold on
loglog(h,H1_err_ts,'c^--','linewidth',1.8, 'MarkerSize',12);
hold on
loglog(h,H1_err_ts_lin,'ks-.','linewidth',1.8, 'MarkerSize',7);
hold on

x=zeros();
y=x;
for i=1:100
    x(i)=0.01*i;
    y(i)=x(i)/7;
end
loglog(x,y,'b-.','linewidth',1.8);
hold on

xlabel({'$h$'},'Interpreter','latex');
ylabel('error');
legend({'$\| u_{ref} - u_{h,h}\|_{1}$',...
        '$\| u_{ref} - u_{H,H}^{h}\|_{1}$',...
        '$\| u_{ref} - \widetilde{u}_{H,H}^{h}\|_{1}$ ','slope=1'},'Interpreter','latex');
set(gca,'XLim',[0.12,max(h)+0.02]);
set(gca,'YLim',[7*10^(-3),4.5*10^(-2)]);
set(gca,'YTick',10^(-2));
set(gca,'FontSize',25);

figure(3) % The running time
dh = (xr-xl)./h;
plot(dh,time_fem(:,1),'mo-','Linewidth',1.8, 'MarkerSize',15);
hold on
plot(dh,time_ts(:,1),'c^--','Linewidth',1.8, 'MarkerSize',12);
hold on
plot(dh,time_ts_lin(:,1),'ks-.','Linewidth',1.8, 'MarkerSize',7);
hold on

xlabel({'$L/h$'},'Interpreter','latex');
ylabel('time(s)');
legend('Standard FEM','Algorithm 5.1','Algorithm 5.2');
set(gca,'XLim',[30.0,80.0]);
% set(gca,'YLim',[0.1,6.5]);
set(gca,'FontSize',25);
