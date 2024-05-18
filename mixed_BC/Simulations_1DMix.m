function Simulations_1DMix(in_K)
clc
close all

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')

%% Numerical set up
par.K = in_K;                    % Number of spatial grid cells
par.L = 1;                       % Domain length
par.Nright = 1;
par.Pright = 1;
par.Uleft = 0;
par.Uright = 0;
x = linspace(0,par.L,par.K);
par.x = linspace(0,par.L,par.K); % Discretise spatial domain
par.dx = par.x(2)-par.x(1);      % Cell size
t0 = 0;                          % Initial time
tf = 2000;                        % Final time
tspan = linspace(t0,tf,100);     % Time span
dt = tspan(2)-tspan(1);

%% Initial conditions - eq.(28)
steadystate = [ones(2*par.K,1); zeros(par.K,1)];
f = 1+0.5*exp(-((par.x)./0.2).^2);
% f1 = 2 + 10*sin(0.1*pi.*par.x);
% f2 = 1 + 0.5*par.x;
% f3 = 2 + 10*cos(0.05*pi.*par.x);
n0 = (f.*ones(1,par.K)).';      % cells
% n0 = ones(par.K,1);
% p0 = (f.*ones(1,par.K)).';      % collagen
p0 = ones(par.K,1);
u0 = zeros(par.K,1);            % displacement
% u0 = sin(pi.*x+pi/2);
% u0 = (f2.*ones(1,par.K)).';
y0 = [n0;p0;u0];

%% Solve with ODE15i
%%% Compute consistent yp0 
res = @(y,yp)(norm(mechanochemical(y,yp,par)));
disp(['residuum of steady state = ' ...
    num2str(res(steadystate,0*steadystate), '%15.10e')]);
opt = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7));
nfixed = zeros(par.K,1);
pfixed = zeros(par.K,1);
ufixed = zeros(par.K,1);
[y0,yp0,resnorm] = decic(@(t,y,yp)(mechanochemical(y,yp,par)), t0, ...
    y0, [nfixed; pfixed; ufixed], ...
    [zeros(3*par.K,1)], zeros(3*par.K,1), opt);
disp(['residuum (from decic) of IC = ' num2str(resnorm, '%15.10e')]);
disp(['residuum (from res()) of IC = ' num2str(res(y0,yp0), '%15.10e')]);

%% Solve
tic
% [t,y] = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,y0,yp0);
sol = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,y0,yp0);
tfinal = sol.stats.tfinal;
t = linspace(t0,tfinal,tfinal/dt+1);
% t = t.';
[y,yp] = deval(sol,t);
toc
%%% Save computed solution to file
filename = ['saved_y1D_' num2str(par.K)];
save(filename, 't', 'y', 'par', 'x', 'yp');

%% Plot
% video_on = true; % Record video: YES (true), NO (false)
% video_filename = [filename '.avi'];
% plot_solution(x,y,yp,t,par,video_on,video_filename);
pic_name = ['res' '.png'];
plot_res(x,y,yp,par,tf,pic_name);

end

%% Main function implementing the model
function f = mechanochemical(y,yp,par)
    %%% Parameter values
    eta = 1;      % viscosity
    E = 1;        % elasticity
    s = 1;        % substrate elasticity
    D = 0.01;     % diffusion
    Dp = 1e-3;       % diffusion for collagen
    an = 0.5;       % cell recruitment rate
    dn = 0.5;       % cell decay rate
    m = 0.1;        % collagen production rate
    dp = 0.1;       % collagen decay rate
    tau = 0.5;    % cell traction
    a1 = 0.5;       % stress related recruitment rate
    
    %%% Reshape input vectors
    [n,p,u] = deal(y(1:par.K),y(par.K+1:2*par.K),y(2*par.K+1:3*par.K));
    ntilde = [n(2); n; par.Nright]; 
    ptilde = [p(2); p; par.Pright]; 
    utilde = [0; u; 0];
    [np,pp,up] = deal(yp(1:par.K),yp(par.K+1:2*par.K),...
         yp(2*par.K+1:3*par.K));
    uptilde = [0; up; 0];

    %%% Traction force term
    % Tr = tau*p.*n;
    n0 = 1.1;
    k2 = 27;
    h1 = (n.^k2)./(n0^k2*ones(size(n)) + n.^k2);
    Tr = tau.*p.*h1;
    % Trtilde = [Tr(2); Tr; tau*ptilde(end)*ntilde(end)];
    % for when traction force term has a hill function
    n2 = (ntilde(end)^k2)/(n0^k2 + ntilde(end)^k2);
    Trtilde = [Tr(2); Tr; tau*ptilde(end)*n2];

    %%% Equation for n
    % Advection velocity at grid cell interfaces - eq.(S.6)
    sig = eta*Mx(uptilde, par) + E*Mx(utilde,par) + Tr;
    % disp(sig.');
    sig0 = 0.2*ones(size(sig));
    k1 = 5;
    fsig = (sig.^k1)./(sig0.^k1 + sig.^k1);
    % fn(n,n',p,u') = 0 - eq.(S.5)
    fn = np - D*Mxx(ntilde, par) + 0*MA1(ntilde, up, par) + MA2(n,up,par) - a1*fsig - an*ones(size(n)) + dn*n;
    % fn = np;

    %%% Equation for p
    % fp(p,p',u') = 0 - eq.(S.10)
    fp = pp - Dp*Mxx(ptilde,par) + 0*MA1(ptilde, up, par) + MA2(p,up,par) - m*n + dp*p;
    % fp = pp;

    %%% Equation for u 
    % Traction term - eq.(S.12)-(S.14)
    % fu(n,n',p,p',u,u') = 0 - eq.(S.11)
    fu = eta*Mxx(uptilde, par) + E*Mxx(utilde,par) + Mx(Trtilde, par) - s*p.*u;
    % fu = up - 0.1*sin(pi.*par.x+pi/2).';
    % fu = up + 0.1;

    %%% Full system - eq.(S.1)
    f = [fn; fp; fu];
end

%% Annexed functions

%%% Compute first order derivative on grid cell interfaces 
%%% second order approximation, central
function dx1 = Mx(y,par)
    persistent Mdx;     % Mdx is a K x K+2 matrix
    c1 = [-1; zeros(par.K-1,1)];
    cn = [zeros(par.K-1,1); 1];
    Mdx = (0.5/par.dx)*[c1, diag(ones(par.K-1,1),1)-diag(ones(par.K-1,1),-1) ,cn];
    dx1 = Mdx*y;
end

%%% Compute second order derivative on grid cell interfaces 
%%% second order approximation, central
function dx2 = Mxx(y,par)
    persistent Mdxx;     % Mdx is a K x K+2 matrix
    c1 = [1; zeros(par.K-1,1)];
    cn = [zeros(par.K-1,1); 1];
    Mdxx = (1/par.dx/par.dx)*[c1, -2*eye(par.K) + diag(ones(par.K-1,1),1) + diag(ones(par.K-1,1),-1) ,cn];
    dx2 = Mdxx*y;
end

%%% Computing second order derivative for displacement, without enforcing
%%% BC on the right
function dx2 = MxxfreeRrightBC(y,par)
    persistent Mdxx;     % Mdx is a K x K+2 matrix
    c1 = [1; zeros(par.K-1,1)];
    cn = [zeros(par.K-1,1); 1];
    Mdxx = (1/par.dx/par.dx)*[c1, -2*eye(par.K) + diag(ones(par.K-1,1),1) + diag(ones(par.K-1,1),-1) ,cn];
    dx2 = Mdxx*y;
    uxx = MxfreeRightBC([0; MxfreeRightBC(y,par); 0],par);
    dx2(end) = uxx(end);
end

% function dx2 = MxxfreeLeftBC(y,par)
%     persistent Mdxx;     % Mdx is a K x K+2 matrix
%     c1 = [1; zeros(par.K-1,1)];
%     cn = [zeros(par.K-1,1); 1];
%     Mdxx = (1/par.dx/par.dx)*[c1, -2*eye(par.K) + diag(ones(par.K-1,1),1) + diag(ones(par.K-1,1),-1) ,cn];
%     dx2 = Mdxx*y;
%     uxx = Mxu([0; Mxu(y,par); 0],par);
%     dx2(end) = uxx(end);
% end

function dx1 = MxfreeRightBC(y,par)
    persistent Mdx;     % Mdx is a K x K+2 matrix
    c1 = [-1; zeros(par.K-1,1)];
    cn = [zeros(par.K-1,1); 1];
    Mdx = (0.5/par.dx)*[c1, diag(ones(par.K-1,1),1)-diag(ones(par.K-1,1),-1) ,cn];
    dx1 = Mdx*y;
    % dx1(1) = (1/par.dx)*(y(3)-y(2));
    dx1(end) = (1/par.dx)*(y(end-1)-y(end-2));
end

%%% Compute advection at grid cell interfaces using first order upwinding
%%% with advective velocity given at grid cell interfaces - def.(S.7)
%%% the current scheme messes up the simplest test of advection if we
%%% enforce the flux at the left boundary. If not, it runs fine until it
%%% reaches the left boundary for a while. 
function fluxdiffx1 = MA1(y, vel, par)
    % compute flux accross cell interfaces using first order upwinding
    rightBC = y(end);
    y = y(2:end-1);
    % flux at K grid interfaces
    flux = NaN(size(vel));

    % take the cell center values for material, using the K by K-1 matrix
    % get K-1 grid cell centers
    c1 = [1; zeros(par.K-2,1)];
    yavg = 0.5*[c1, eye(par.K-1)+diag(ones(par.K-2,1),-1)]*y;

    % interate over all but boundary cell interfaces
    for i = 2:size(flux)-1
        if (vel(i)>0)
            flux(i) = vel(i)*yavg(i-1);
        else 
            flux(i) = vel(i)*yavg(i);
        end
    end

    % for flux at the left Neumann boundary
    if vel(1)<0
        flux(1) = vel(1)*yavg(1);
    else
        flux(1) = vel(1)*y(1);
    end

    % for flux at the right boundary
    if vel(end)>0
        flux(end) = vel(end)*yavg(end);
    else
        flux(end) = vel(end)*rightBC;
    end

    % compute flux difference per grid cell - def.(S.7) and (S.4)
    fluxdiffx1 = Mx([0; flux; 0],par);
    fluxdiffx1(1) = (flux(2)-flux(1))/par.dx;
    fluxdiffx1(end) = (flux(end)-flux(end-1))/par.dx;
    % fluxdiffx1(end) = 0;
    
    % compute flux difference per grid cell via first order FD, assuming
    % flux gradient at the right boundary is 0
    % fluxdiffx1 = [(1/par.dx)*(flux(2:end)-flux(1:(end-1))); 0];

    % compute flux difference per grid cell via first order FD, assuming
    % flux gradient at the right boundary is the same as to its left 
    % fluxdiffx1 = (1/par.dx)*(flux(2:end)-flux(1:(end-1)));
    % fluxdiffx1 = [fluxdiffx1; fluxdiffx1(end)];
end

%%% Compute flux gradient using first order upwinding without center
%%% averaging
%%% Take material at K grid cells
function fluxdiffx1 = MA2(y, vel, par)
    flux = zeros(size(y));

    for i=2:size(flux)
        if vel(i)>0
            flux(i) = y(i-1)*vel(i);
        else
            flux(i) = y(i)*vel(i);
        end
    end

    % flux(1) = 0;

    flux(1) = vel(1)*y(1);

    % compute flux difference per grid cell - def.(S.7) and (S.4)
    % fluxdiffx1 = Mx([0; flux; 0],par);
    % fluxdiffx1(1) = (flux(2)-flux(1))/par.dx;
    % fluxdiffx1(end) = (flux(end)-flux(end-1))/par.dx;
    
    % compute flux difference per grid cell via first order FD, assuming
    % flux gradient at the right boundary is 0
    fluxdiffx1 = [(1/par.dx)*(flux(2:end)-flux(1:(end-1))); 0];

    % compute flux difference per grid cell via first order FD, assuming
    % flux gradient at the right boundary is the same as to its left 
    % fluxdiffx1 = (1/par.dx)*(flux(2:end)-flux(1:(end-1)));
    % fluxdiffx1 = [fluxdiffx1; fluxdiffx1(end)];
end

% Take cells or collagen without ghost points (K cell interfaces)
% function fluxdiffx1 = MA2(y, vel, par)
%     flux = zeros(size(y));
%     % velocity shifted to grid centers
%     c1 = [1; zeros(par.K-2,1)];
%     vel = 0.5*[c1, eye(par.K-1)+diag(ones(par.K-2,1),-1)]*vel;
% 
%     % compute flux at the grid cell interfaces
%     for i = 1:size(vel)
%         if vel(i)>0
%             flux(i) = flux(i)-y(i)*vel(i);
%             flux(i+1) = flux(i+1)+y(i)*vel(i);
%         else
%             flux(i) = flux(i)-y(i+1)*vel(i);
%             flux(i+1) = flux(i+1)+y(i+1)*vel(i);
%         end
%     end
% 
%     % compute flux difference per grid cell - def.(S.7) and (S.4)
%     fluxdiffx1 = Mx([flux(2); flux; 0], par);
%     fluxdiffx1(1) = (flux(2)-flux(1))/par.dx;
%     fluxdiffx1(end) = (flux(end)-flux(end-1))/par.dx;
% end

%%% Plot solution 
function plot_solution(x,y,yp,t,par,video_on,video_filename)
    if video_on % Initialise video
        vid = VideoWriter(video_filename);
        open(vid);
        figure('Units','normalized','Position',[0 0 0.5 0.45])
    end
    maxu = 10^(-6);

    for i=1:length(t)
        clf
        n = [y(1:par.K,i)];
        p = [y(par.K+1:2*par.K,i)];
        u = [y(2*par.K+1:3*par.K,i)];
        v = [yp(2*par.K+1:3*par.K,i)];
        if max(abs(u))>maxu
          maxu = max(abs(u))+0.1*max(abs(u));
        end
        subplot(1,4,1)
        plot(x,n)
        title('$n(t,x)$')
        % ylim([0.8,max(1.2,max(n))])
        axis square
        subplot(1,4,2)
        plot(x,p)
        title('$\rho(t,x)$')
        % ylim([0.8,max(1.2,max(p))])
        axis square
        subplot(1,4,3)
        plot(x,u)
        title('$u(t,x)$')
        axis square
        subplot(1,4,4)
        plot(x,v)
        title('$v(t,x)$')
        axis square
        a = axes;
        t1 = title([' (t=',num2str(t(i)),')']);
        a.Visible = 'off'; 
        t1.Visible = 'on'; 
        drawnow
        if video_on % Record video
            frame = getframe(gcf);
            size(frame.cdata);
            writeVideo(vid,frame);
            pause(0.1)
        end
    end
    if video_on % Close video
        close(vid)
    end
end

%%% Plot solution 
function plot_res(x,y,yp,par,tf,pic_name)
    figure('Units','normalized','Position',[0 0 0.5 0.45])
    n = [y(1:par.K,end)];
    p = [y(par.K+1:2*par.K,end)];
    u = [y(2*par.K+1:3*par.K,end)];
    v = [yp(2*par.K+1:3*par.K,end)];
    % res_file = 'res';
    % save(res_file,"n","p","u");
    subplot(1,4,1)
    plot(x,n)
    title('$n(t,x)$')
    % ylim([0.8,max(1.2,max(n))])
    axis square
    subplot(1,4,2)
    plot(x,p)
    title('$\rho(t,x)$')
    % ylim([0.8,max(1.2,max(p))])
    axis square
    subplot(1,4,3)
    plot(x,u)
    title('$u(t,x)$')
    axis square
    subplot(1,4,4)
    plot(x,v)
    title('$v(t,x)$')
    axis square
    a = axes;
    a.Visible = 'off';
    t1 = title([' (t=',num2str(tf),')']);
    t1.Visible = "on";
    saveas(gcf,pic_name);
end